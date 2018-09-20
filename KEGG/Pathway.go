package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
)

const HELP = `
KEGG pathway process, usage:

1. update local data table ($EXECUTINGPATH/KEGG_data/KEGG_organism.tsv):
    $ Pathway  Update

2. download organisms keg file(s):
    $ Pathway  Get  hsa mmu ath

3. get keg file of an organism from local:
    $ Pathway  get  hsa
    Note: make sure you have download organisms' keg files and achieve to 
    $EXECUTINGPATH/KEGG_data/Pathway_keg.tar

4. find match species name or code in local data table:
    $ Pathway  match  "Rhinopithecus roxellana"
    $ Pathway  match  Rhinopithecus+roxellana
    $ Pathway  match  rro

5. download pathway html:
    $ Pathway  HTML  hsa00001.keg.gz  ./hsa00001
    Note: existing html files will not be overwritten

6. convert keg format to tsv (file or stdout):
    $ Pathway  tsv  hsa00001.keg.gz  hsa00001.keg.tsv

    output tsv header:
    C_id  id  gene  KO  EC

7. download species keg, convert to tsv and download html files:
    $ Pathway  species  Rhinopithecus+roxellana
    Note: existing html files will be overwritten

author: d2jvkpn
version: 0.8.5
release: 2018-09-20
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main() {
	nargs := len(os.Args) - 1

	if nargs == 0 || os.Args[1] == "-h" || os.Args[1] == "--help" {
		fmt.Println(HELP)
		return
	}

	cmd := os.Args[1]
	ep, _ := exec.LookPath(os.Args[0])
	datatsv := filepath.Dir(ep) + "/KEGG_data/KEGG_organism.tsv"

	switch {
	case cmd == "Update" && nargs == 1:
		Update(datatsv)

	case cmd == "Get" && nargs > 1:
		Get(os.Args[2:])

	case cmd == "get" && nargs == 2:
		ok := Get_local(os.Args[2] + "00001.keg.gz", 
			filepath.Dir(ep)+"/KEGG_data/Pathway_keg.tar")

		if !ok {
			os.Exit(1)
		}

	case cmd == "HTML" && nargs == 3:
		DownloadHTML(os.Args[2], os.Args[3], false)

	case cmd == "match" && nargs == 2:
		record, found := Match(os.Args[2], datatsv)

		if found {
			fmt.Printf("Entry: %s\nCode: %s\nSpecies: %s\nLineage: %s\n",
				record[0], record[1], record[2], record[3])
		} else {
			fmt.Fprintln(os.Stderr, "NotFound")
		}

	case cmd == "tsv" && (nargs == 3 || nargs == 2):
		if nargs == 3 {
			ToTSV(os.Args[2], os.Args[3])
		} else {
			ToTSV(os.Args[2], "")
		}

	case cmd == "species" && nargs == 2:
		record, found := Match(formatSpeciesName(os.Args[2]), datatsv)

		if found {
			fmt.Printf("Entry: %s\nCode: %s\nSpecies: %s\nLineage: %s\n",
				record[0], record[1], record[2], record[3])

			log.Printf("Querying %s\n", record[1]+"00001.keg")

			ok := getkeg(record[1] + "00001.keg")
			if !ok {
				os.Exit(1)
			}

			ToTSV(record[1]+"00001.keg.gz", record[1]+"00001.keg.tsv")

			DownloadHTML(record[1]+"00001.keg.gz", record[1]+"00001", true)
		} else {
			fmt.Fprintln(os.Stderr, "NotFound")
		}

	default:
		log.Fatal(HELP)
	}
}

func formatSpeciesName(name string) string {
	wds := strings.Fields(strings.Replace(name, "+", " ", -1))
	re := regexp.MustCompile("[A-Za-z][a-z]+")

	for i, _ := range wds {
		if re.Match([]byte(wds[i])) {
			wds[i] = strings.ToLower(wds[i])
		}
	}

	wds[0] = strings.Title(wds[0])
	return (strings.Join(wds, " "))
}

func DownloadHTML(keg, outdir string, overwrite bool) {
	var bts []byte
	var err error

	frd, err := os.Open(keg)
	if err != nil {
		log.Fatal(err)
	}
	defer frd.Close()

	err = os.MkdirAll(outdir, 0755)
	if err != nil {
		log.Fatal(err)
	}

	log.Printf("Save html files to %s\n", outdir)

	if strings.HasSuffix(keg, ".gz") {
		var gz *gzip.Reader
		gz, err = gzip.NewReader(frd)
		if err != nil {
			return
		}

		bts, err = ioutil.ReadAll(gz)
	} else {
		bts, err = ioutil.ReadAll(frd)
	}

	re := regexp.MustCompile("PATH:[a-z]+[0-9]+")
	PATHs := re.FindAllString(string(bts), -1)

	for i, _ := range PATHs {
		PATHs[i] = strings.TrimLeft(PATHs[i], "PATH:")
	}

	var wg sync.WaitGroup
	ch := make(chan struct{}, 10)

	for _, p := range PATHs {
		ch <- struct{}{}
		wg.Add(1)
		go gethtml(p, outdir, overwrite, ch, &wg)
	}
	wg.Wait()
}

func gethtml(p, outdir string, overwrite bool, ch <-chan struct{},
	wg *sync.WaitGroup) {

	defer func() { <-ch }()
	defer wg.Done()

	html := outdir + "/" + p + ".html"
	png := outdir + "/" + p + ".png"
	url := "http://www.genome.jp/kegg"
	code := p[:(len(p) - 5)]

	if _, err := os.Stat(html); err == nil && !overwrite {
		return
	}

	htmlurl := fmt.Sprintf(url+"-bin/show_pathway?%s", p)
	pngurl := fmt.Sprintf(url+"/pathway/%s/%s.png", code, p)

	htmlresp, err := http.Get(htmlurl)
	if err != nil {
		log.Println(err)
		return
	}
	defer htmlresp.Body.Close()

	htmlbody, err := ioutil.ReadAll(htmlresp.Body)
	if err != nil {
		log.Println(err)
		return
	}

	text := string(htmlbody)
	if !strings.HasSuffix(text, "</html>\n") {
		return
	}

	re, _ := regexp.Compile("\\[[\\S\\s]+?\\]")
	text = re.ReplaceAllString(text, "")

	re, _ = regexp.Compile("\\<script[\\S\\s]+?\\</script\\>")
	text = re.ReplaceAllString(text, "")

	re, _ = regexp.Compile("\\<style[\\S\\s]+?\\</style\\>")
	text = re.ReplaceAllString(text, "")

	re, _ = regexp.Compile("\\<link[\\S\\s]+?\\/\\>[ \n]+")
	text = re.ReplaceAllString(text, "\n")

	re, _ = regexp.Compile("\\<table[\\S\\s]+?\\</table\\>")
	text = re.ReplaceAllString(text, "")

	re, _ = regexp.Compile("\\<div[\\S\\s]+?\\</form>")
	text = re.ReplaceAllString(text, "</div>")

	re, _ = regexp.Compile("\\<body\\>[ \n]+")
	text = re.ReplaceAllString(text, "<body>\n<div align=\"center\">\n")

	text = strings.Replace(text, fmt.Sprintf("/kegg/pathway/%s/", code), "", 1)

	text = strings.Replace(text, "/dbget-bin/www_bget?",
		"https://www.genome.jp/dbget-bin/www_bget?", -1)

	pngresp, err := http.Get(pngurl)
	if err != nil {
		log.Println(err)
		return
	}
	defer pngresp.Body.Close()

	pngbody, err := ioutil.ReadAll(pngresp.Body)
	if err != nil {
		log.Println(err)
		return
	}

	err = ioutil.WriteFile(png, pngbody, 0664)
	if err != nil {
		log.Println(err)
		return
	}

	err = ioutil.WriteFile(html, []byte(text), 0664)
	if err != nil {
		log.Println(err)
		return
	}
}

func ToTSV(keg, tsv string) {
	scanner, frd, err := ReadInput(keg)
	if err != nil {
		log.Fatal(err)
	}
	defer frd.Close()

	type Write interface {
		Write(p []byte) (n int, err error)
	}

	var TSV Write

	if tsv != "" {
		err = os.MkdirAll(filepath.Dir(tsv), 0755)
		if err != nil {
			log.Fatal(err)
		}

		fwt, err := os.Create(tsv)
		if err != nil {
			log.Fatal(err)
		}

		defer fwt.Close()
		TSV = fwt
	} else {
		TSV = os.Stdout
	}

	var line string
	var fds [5]string

	TSV.Write([]byte("C_id\tid\tgene\tKO\tEC\n"))
	A, B := "", ""

	for scanner.Scan() {
		line = scanner.Text()
		if len(line) < 2 {
			continue
		}

		switch line[0] {
		case 'A':
			A = strings.Replace(line, " ", ":", 1)

		case 'B':
			B = strings.Replace (strings.Replace(line, "  ", "", 1),
				" ", ":", 1)

		case 'C':
			fds[3], fds[4] = A, B

			tmp := strings.SplitN(strings.TrimLeft(line, "C    "), " ", 2)
			fds[0], fds[1] = "C" + tmp[0], ""
			fds[2] = strings.TrimRight (tmp[1], "]")

			if strings.Contains(fds[2], " [") {
				copy(fds[1:3], strings.SplitN(fds[2], " [", 2))
				fds[1], fds[2] = fds[2], fds[1]
			}

			TSV.Write([]byte("#" + strings.Join (fds[0:], "\t") + "\n"))

		case 'D':
			sep := "\t"

			if !strings.Contains(line, "\t") {
				sep = "; " // for KAAS output keg
			}

			tmp := strings.SplitN(strings.TrimLeft(line, "D      "), sep, 2)
			if len(tmp) != 2 {
				continue
			}

			if strings.Contains(tmp[0], " ") {
				copy(fds[1:3], strings.SplitN(tmp[0], " ", 2))
			} else {
				fds[1], fds[2] = tmp[0], "" // for KAAS output keg
			}

			tmp[1] = strings.Replace(tmp[1], "  ", " ", 1) // for KAAS output keg

			if strings.Contains(tmp[1], " [EC:") {
				copy(fds[3:5], strings.SplitN(
					strings.Replace(tmp[1], " [EC:", "\t[EC:", 1), "\t", 2))
			} else {
				fds[3], fds[4] = tmp[1], ""
			}

			TSV.Write([]byte(strings.Join (fds[0:], "\t") + "\n"))

		default:
			continue

		}
	}

	if tsv != "" {
		log.Printf("Converted %s to %s\n", keg, tsv)
	}
}

func Match(name, datatsv string) (record []string, ok bool) {
	file, err := os.Open(datatsv)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	species := formatSpeciesName(name)
	scanner := bufio.NewScanner(file)
	scanner.Scan() // skip header

	for scanner.Scan() {
		record = strings.Split(scanner.Text(), "\t")

		ok = (name == record[1] ||
			species == strings.Split(record[2], " (")[0])

		if ok {
			return
		}
	}

	return
}

func Update(saveto string) {
	log.Printf("Update %s\n", saveto)

	resp, err := http.Get("http://rest.kegg.jp/list/organism")
	if err != nil {
		log.Println(err)
		return
	}
	defer resp.Body.Close()

	body, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		log.Println(err)
		return
	}

	err = os.MkdirAll(filepath.Dir(saveto), 0755)
	if err != nil {
		log.Fatal(err)
	}

	fwt, err := os.Create(saveto)
	if err != nil {
		log.Println(err)
		return
	}
	defer fwt.Close()

	fwt.Write([]byte("Entry\tCode\tSpecies\tLineage\n"))
	fwt.Write(body)
}

func Get_local(keg, path string) (ok bool) {
	Cmd := exec.Command("tar", "-xf", path, keg)
	err := Cmd.Run()

	if err != nil {
		log.Printf("Failed to get %s from %s\n", keg, path)
	} else {
		ok = true
	}

	return ok
}

func Get(codes []string) {
	ch := make(chan struct{}, 10)
	var wg sync.WaitGroup

	log.Printf("Request organism code(s): %s\n", strings.Join(codes, " "))

	for _, v := range codes {
		ch <- struct{}{}
		wg.Add(1)
		go func(p string, ch <-chan struct{}, wg *sync.WaitGroup) {
			defer func() { <-ch }()
			defer wg.Done()
			getkeg(p)
		}(v+"00001.keg", ch, &wg)
	}

	wg.Wait()
}

func getkeg(p string) (ok bool) {
	resp, err := http.Get(fmt.Sprintf("http://www.kegg.jp/kegg-bin"+
		"/download_htext?htext=%s&format=htext&filedir=", p))

	if err != nil {
		log.Println(err)
		return
	}
	defer resp.Body.Close()

	body, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		log.Println(err)
		return
	}

	lines := strings.Split(string(body), "\n")

	if !strings.HasPrefix(lines[len(lines)-2], "#Last updated:") {
		log.Printf("Failed to get %s\n", p)
		return
	}

	file, err := os.Create(p + ".gz")
	if err != nil {
		log.Println(err)
		return
	}
	defer file.Close()

	gw := gzip.NewWriter(file)
	gw.Write(body)
	gw.Close()
	log.Printf("Saved %s.gz\n", p)

	ok = true
	return
}

func ReadInput(s string) (scanner *bufio.Scanner, file *os.File, err error) {
	file, err = os.Open(s)
	if err != nil {
		return
	}

	if strings.HasSuffix(s, ".gz") {
		var gz *gzip.Reader
		gz, err = gzip.NewReader(file)
		if err != nil {
			return
		}

		scanner = bufio.NewScanner(gz)
	} else {
		scanner = bufio.NewScanner(file)
	}

	return
}
