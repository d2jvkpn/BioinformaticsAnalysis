package main

import (
	"os"
	"fmt"
	"log"
	"io/ioutil"
	"strings"
	"sync"
	"bufio"
	"compress/gzip"
	"net/http"
	"path/filepath"
	"regexp"
)

const HELP = `
Update local data table for command "match":
    $ Pathway  update

Find match species in local data table:
    $ Pathway  match  "Rhinopithecus roxellana"
    $ Pathway  match  Rhinopithecus+roxellana

Get organisms keg file:
    $ Pathway  get  hsa mmu ath

Convert keg format to tsv:
    $ Pathway  totsv  hsa00001.keg.gz  hsa00001.keg.tsv
    $ Pathway  totsv  hsa00001.keg  hsa00001.keg.tsv

Get species keg and convert to tsv:
    $ Pathway  species  Rhinopithecus+roxellana

author: d2jvkpn
version: 0.2
release: 2018-09-08
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

const URL = "http://www.kegg.jp/kegg-bin/download_htext?htext=%s&format=htext&filedir="

func main () {
	nargs := len (os.Args) - 1

	if nargs == 0 || os.Args[1] == "-h" ||  os.Args[1] == "--help" {
		fmt.Println (HELP); return
	}

	cmd := os.Args[1]
	datatsv := filepath.Dir (os.Args[0]) + "/data/KEGG_organism.tsv"

	switch {
	case nargs == 1 && cmd == "update":
		Update (datatsv)

	case nargs > 1 && cmd == "get":
		Get (os.Args[2:])

	case nargs == 2 && cmd == "match":
		record, found := Match (formatSpeciesName (os.Args[2]), datatsv)

		if found {
			fmt.Printf ("Entry: %s\nCode: %s\nSpecies: %s\nLineage: %s\n\n",
			record[0], record[1], record[2], record[3])
		} else {
			fmt.Println ("NotFound")
		}

	case nargs == 3 && cmd == "totsv":
		ToTSV (os.Args[2], os.Args[3])

	case nargs == 2 && cmd == "species":
		record, found := Match (formatSpeciesName (os.Args[2]), datatsv)

		if found {
			Get ( []string {record[1]} )
			ToTSV (record[1] + "00001.keg.gz", record[1] + "00001.keg.tsv")
		} else {
			log.Fatal ("NotFound")
		}

	default:
		log.Fatal (HELP)
	}
}


func formatSpeciesName (name string) string {
	wds := strings.Fields (strings.Replace (name, "+", " ", -1))
	re := regexp.MustCompile("[A-Za-z][a-z]+")

	for i, _ := range wds {
		if re.Match ([]byte(wds[i])) { wds[i] = strings.ToLower(wds[i]) }
	}
	
	wds[0] = strings.Title (wds[0])
	return (strings.Join (wds, " "))
}


func ToTSV (keg, tsv string) {
	scanner, frd, err := ReadInput (keg)
	if err != nil { log.Fatal (err) }
	defer frd.Close()

	TSV, err := os.Create (tsv)
	if err != nil { log.Fatal (err) }
	defer TSV.Close()

	var line string
	var fds [11]string

	TSV.Write ([] byte ("C_id\tC_entry\tC_name\tgene_id\tgene\t" + 
	"A_id\tA_name\tB_id\tB_name\tKO\tEC\n"))

	for scanner.Scan () {
		line = scanner.Text()
		if len (line) < 2 { continue }

		switch line[0] {
		case 'A':
			copy (fds[5:7], strings.SplitN (line, " ", 2))

		case 'B':
			copy (fds[7:9], strings.SplitN (strings.TrimLeft (line, "B  "),
			" ", 2) )

			fds[7] = "B" + fds[7]

		case 'C':
			tmp := strings.SplitN (strings.TrimLeft (line, "C	"), " ", 2)
			fds[0], fds[1] = "C" + tmp[0], ""
			fds[2] = strings.TrimRight (tmp[1], "]") 

			if strings.Contains (fds[2], " [") {
				copy (fds[1:3], strings.SplitN (fds[2], " [", 2))
				fds[1], fds[2] = fds[2], fds[1]
			}

			fds[3], fds[4], fds[9], fds[10] = "", "", "", ""

			TSV.Write ([]byte (strings.Join (fds[0:], "\t") + "\n"))

		case 'D':
			tmp := strings.SplitN (strings.TrimLeft (line, "D	  "), "\t", 2)
			if len(tmp) != 2 { continue }

			copy (fds[3:5], strings.SplitN (tmp[0], " ", 2))
			copy (fds[9:11], strings.SplitN (tmp[1], " ", 2))

			TSV.Write ([]byte (strings.Join ( []string {fds[0], "-", "-", 
			fds[3], fds[4], "-", "-", "-", "-", fds[9], fds[10] }, "\t") + "\n"))

		default:
			continue

		}
	}
}


func Match (name, datatsv string) (record []string, ok bool) {
	file, err := os.Open (datatsv)
	if err != nil { log.Fatal(err) }
	defer file.Close()

	scanner := bufio.NewScanner (file)
	scanner.Scan () // skip header

	for scanner.Scan () {
		record = strings.Split (scanner.Text(), "\t")
		ok = (name == strings.Split (record[2], " (")[0])
		if ok { return }
	}

	return
}


func Update (saveto string) {
	resp, err := http.Get ("http://rest.kegg.jp/list/organism")
	if err != nil { log.Println (err); return }
	defer resp.Body.Close ()

	body, err := ioutil.ReadAll (resp.Body)
	if err != nil { log.Println (err); return }

	err = os.MkdirAll (filepath.Dir (saveto), 0755)
	if err != nil { log.Fatal (err)}

	fwt, err := os.Create (saveto)
	if err != nil { log.Println (err); return }
	defer fwt.Close ()

	fwt.Write ([]byte ("Entry\tCode\tSpecies\tLineage\n"))
	fwt.Write (body)
}


func Get (codes []string) {
	ch := make (chan struct {}, 10)
	var wg sync.WaitGroup

	log.Printf ("Request organism code(s):\n	%s\n", 
	strings.Join (codes, " "))

	for _, v := range codes {
		ch <- struct{}{}
		wg.Add (1)
		go getkeg (v + "00001.keg", ch, &wg)
	}

	wg.Wait ()
} 


func getkeg (p string, ch <- chan struct{}, wg *sync.WaitGroup) {
	defer func () { <- ch }()
	defer wg.Done ()

	// log.Printf ("Querying %s...\n", p)

	resp, err := http.Get (fmt.Sprintf (URL, p))
	if err != nil { log.Println (err); return }
	defer resp.Body.Close ()

	body, err := ioutil.ReadAll (resp.Body)
	if err != nil { log.Println (err); return }

	lines := strings.Split (string (body), "\n")

	if ! strings.HasPrefix (lines[len (lines)-2], "#Last updated:") {
		log.Printf ("Failed to get %s\n", p)
		return
	}

	file, err := os.Create (p + ".gz")
	if err != nil { log.Println (err); return }
	defer file.Close ()

	gw := gzip.NewWriter (file)
	gw.Write (body)
	gw.Close ()

	log.Printf ("Saved %s...\n", p)
}


func ReadInput (s string) (scanner *bufio.Scanner, file *os.File, err error) {
	file, err = os.Open (s)
	if err != nil { return } 

	if strings.HasSuffix (s, ".gz") {
		var gz *gzip.Reader
		gz, err = gzip.NewReader (file)
		if err != nil { return }

		scanner = bufio.NewScanner (gz)
	} else {
		scanner = bufio.NewScanner (file)
	}
 
	return
}
