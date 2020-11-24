package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/PuerkitoBio/goquery"
	"log"
	"net/http"
	"net/url"
	"os"
	"path"
	"regexp"
	"strings"
	"time"
)

const USAGE = `Query NCBI genomic information by providing taxonomy id, 
scientific name or genome url, usage:
  $ biocrawler_NCBIgenome  [-d ./]  <input>
`

const LISENSE = `
author: d2jvkpn
version: 0.7
release: 2018-09-04
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main() {
	var dir string
	flag.StringVar(&dir, "d", "./", "save result to directory")

	flag.Usage = func() {
		fmt.Println(USAGE)
		flag.PrintDefaults()
		fmt.Println(LISENSE)
		os.Exit(2)
	}

	flag.Parse()

	if flag.NArg() != 1 {
		flag.Usage()
	}

	input := flag.Args()[0]
	prefix := "https://www.ncbi.nlm.nih.gov/genome"

	if yes, _ := regexp.MatchString("^[1-9][0-9]*$", input); yes {
		input = fmt.Sprintf(prefix+"/?term=txid%s[orgn]", input)

	} else if !strings.Contains(input, prefix) {
		input = strings.Join(strings.Fields(input), " ")
		input = prefix + "/?term=" + url.QueryEscape(input)
	}

	log.Println("querying", input)

	var result Genomic
	var err error
	if err = result.Query(input); err != nil {
		log.Fatal(err)
	}

	if dir == "" {
		fmt.Println(result.ToIni())
	} else {
		err = result.Save(dir)
	}

	if err != nil {
		log.Fatal(err)
	}
}

type Genomic struct {
	Meta, Lineage []string
	Reference     []string
}

func (result *Genomic) Query(query string) (err error) {
	var res *http.Response
	var doc *goquery.Document
	var lineage *goquery.Selection

	if res, err = http.Get(query); err != nil {
		return
	}
	defer res.Body.Close()

	if res.StatusCode != 200 {
		err = fmt.Errorf("%d %s", res.StatusCode, res.Status)
		return
	}

	if doc, err = goquery.NewDocumentFromReader(res.Body); err != nil {
		return
	}

	if strings.Contains(doc.Find("#messagearea").Text(), "No items found") {
		err = errors.New("genome not found")
		return
	}

	result.Meta = append(result.Meta, "URL\t"+query)

	result.Meta = append(result.Meta,
		time.Now().Format("AcessTime\t2006-01-02 15:04:05 -0700"))

	result.Meta = append(result.Meta, "Name\t")

	lineage = doc.Find("span.GenomeLineage")

	if lineage == nil {
		lineage = doc.Find("SPAN[CLASS=GenomeLineage]")
	}

	lineage.Eq(0).Find("a").Each(func(i int,
		sel *goquery.Selection) {

		if i%2 == 0 {
			href, _ := sel.Attr("href")
			result.Lineage = append(result.Lineage,
				path.Base(href)+"\t"+sel.Text())
		}
	})

	if len(result.Lineage) == 0 {
		err = errors.New("species lineage information not found")
		return
	}

	sname := strings.Split(result.Lineage[len(result.Lineage)-1], "\t")[1]
	result.Meta[2] = result.Meta[2] + "NCBI__" + url.QueryEscape(sname)

	doc.Find("div.refgenome_sensor").Eq(0).Find("a").Each(func(i int,
		sel *goquery.Selection) {

		href, _ := sel.Attr("href")

		if strings.HasPrefix(href, "/") {
			href = "https://www.ncbi.nlm.nih.gov" + href
		}

		if sel.Text() == "genome" && strings.HasPrefix(href, "ftp:") {
			result.Meta[2] += "__" +
				strings.TrimRight(path.Base(href), "_genomic.fna.gz")
		}

		result.Reference = append(result.Reference, sel.Text()+"\t"+href)
	})

	return
}

func (result *Genomic) Save(dir string) (err error) {
	var pd, ftp, cmd string
	var fds []string
	var wtj, wts *os.File

	pd = path.Join(dir, strings.SplitN(result.Meta[2], "\t", 2)[1])

	if err = os.MkdirAll(pd, 0755); err != nil {
		return
	}

	if wtj, err = os.Create(pd + "/NCBI_genomic.ini"); err != nil {
		return
	}

	defer wtj.Close()

	fmt.Fprintln(wtj, result.ToIni())

	if wts, err = os.Create(pd + "/download.sh"); err != nil {
		return
	}

	defer wts.Close()

	wts.Write([]byte(fmt.Sprintf(
		"#! /bin/bash\n## URL: %s\nset -eu\ncode=0\n",
		strings.SplitN(result.Meta[0], "\t", 2)[1]),
	))

	cmd = `
{ test -f $(dirname $0)/%[1]s.wget.failed && rm $(dirname $0)/%[1]s.wget.failed
wget -c -O $(dirname $0)/%[1]s -o $(dirname $0)/%[1]s.wget.logging \
%[2]s \
&& gzip -t $(dirname $0)/%[1]s && rm $(dirname $0)/%[1]s.wget.logging ||
{ mv $(dirname $0)/%[1]s.wget.logging $(dirname $0)/%[1]s.wget.failed; code=1; }
} &
`

	for i, _ := range result.Reference {
		ftp = strings.SplitN(result.Reference[i], "\t", 2)[1]

		if !strings.HasPrefix(ftp, "ftp://") ||
			strings.HasSuffix(ftp, "/") {
			continue
		}

		fds = strings.SplitN(path.Base(ftp), "_", -1)

		wts.Write([]byte(fmt.Sprintf(cmd, fds[len(fds)-1], ftp)))
	}

	wts.Write([]byte("\nwait\nexit $code\n"))
	log.Printf("saved %s/{NCBI_genomic.ini,download.sh}\n", pd)

	return
}

func (result *Genomic) ToIni() (str string) {
	str = "\n[Meta]\n"

	for _, s := range result.Meta {
		str += strings.Replace(s, "\t", " = ", 1) + "\n"
	}

	str += "\n[Lineage]\n"
	for _, s := range result.Lineage {
		str += strings.Replace(s, "\t", " = ", 1) + "\n"
	}

	str += "\n[Reference]\n"
	for _, s := range result.Reference {
		str += strings.Replace(s, "\t", " = ", 1) + "\n"
	}

	return
}
