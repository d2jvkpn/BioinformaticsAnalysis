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
	"os/exec"
)

const HELP = `
Update local data table for command "match":
    $ Pathway  update

Find match species name or code in local data table:
    $ Pathway  match  "Rhinopithecus roxellana"
    $ Pathway  match  Rhinopithecus+roxellana
    $ Pathway  match  rro

Get organisms a single keg file from local:
    $ Pathway  get  hsa 

Get organisms keg file(s) by download:
    $ Pathway  Get  hsa mmu ath

Convert keg format to tsv:
    $ Pathway  totsv  hsa00001.keg.gz  hsa00001.keg.tsv

    hsa00001.keg.tsv tsv header:
    C_id C_entry C_name gene_id gene A_id A_name B_id B_name KO KO_name EC

    other output files: gene_id2pathway.tsv pathway.infor.tsv

Get species keg (from local) and convert to tsv:
    $ Pathway  species  Rhinopithecus+roxellana
    output files: rro00001.keg.gz rro00001.keg.tsv

author: d2jvkpn
version: 0.5
release: 2018-09-11
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
	ep, _ := exec.LookPath (os.Args[0])
	datatsv := filepath.Dir (ep) + "/KEGG_data/KEGG_organism.tsv"

	switch {
	case nargs == 1 && cmd == "update":
		Update (datatsv)

	case nargs > 1 && cmd == "get":
		ok := Get_local (os.Args[2], filepath.Dir (ep) + "/KEGG_data/Pathway_kegs")
		if !ok { os.Exit(1) }

	case nargs > 1 && cmd == "Get":
		Get (os.Args[2:])

	case nargs == 2 && cmd == "match":
		record, found := Match (os.Args[2], datatsv)

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
			fmt.Printf ("Entry: %s\nCode: %s\nSpecies: %s\nLineage: %s\n\n",
			record[0], record[1], record[2], record[3])

			ok := Get_local (record[1], 
			filepath.Dir (ep) + "/KEGG_data/Pathway_kegs")

			if !ok { os.Exit(1) }
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
	log.Printf ("Convert %s to tsv\n", keg)

	scanner, frd, err := ReadInput (keg)
	if err != nil { log.Fatal (err) }
	defer frd.Close()

	err = os.MkdirAll (filepath.Dir (tsv), 0755)
	if err != nil { log.Fatal (err)}

	TSV, err := os.Create (tsv)
	if err != nil { log.Fatal (err) }
	defer TSV.Close()

	G2P, err := os.Create (filepath.Dir (tsv) + "/gene_id2pathway.tsv")
	if err != nil { log.Fatal (err) }
	defer G2P.Close()

	Pinfor, err := os.Create (filepath.Dir (tsv) + "/pathway.infor.tsv")
	if err != nil { log.Fatal (err) }
	defer Pinfor.Close()

	var line string
	var fds [12]string
	var isPathway bool

	TSV.Write ([] byte ("C_id\tC_entry\tC_name\tgene_id\tgene\t" + 
	"A_id\tA_name\tB_id\tB_name\tKO\tKO_name\tEC\n"))

	G2P.Write ([]byte ("gene_id\tpathway\tEC\tgene_name\n"))
	Pinfor.Write ([]byte ("pathway\tpathway_name\tA\tB\n"))

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
			tmp := strings.SplitN (strings.TrimLeft (line, "C    "), " ", 2)
			fds[0], fds[1] = "C" + tmp[0], ""
			fds[2] = strings.TrimRight (tmp[1], "]") 

			if strings.Contains (fds[2], " [") {
				copy (fds[1:3], strings.SplitN (fds[2], " [", 2))
				fds[1], fds[2] = fds[2], fds[1]
			}

			fds[3], fds[4], fds[9], fds[10], fds[11] = "", "", "", "", ""
			TSV.Write ([]byte (strings.Join (fds[0:], "\t") + "\n"))
			isPathway = strings.HasPrefix (fds[1], "PATH:")

			if isPathway {
				Pinfor.Write ([]byte ( strings.TrimLeft (fds[1], "PATH:") +
				"\t" + fds[2] + "\t" + fds[6] + "\t" + fds[8] + "\n"))
			}

		case 'D':
			tmp := strings.SplitN (strings.TrimLeft (line, "D	  "), "\t", 2)
			if len(tmp) != 2 { continue }

			copy (fds[3:5], strings.SplitN (tmp[0], " ", 2))
			copy (fds[9:11], strings.SplitN (tmp[1], " ", 2))

			if strings.Contains (fds[10], " [EC:") {
				copy (fds[10:12], strings.SplitN (
				strings.Replace (fds[10], " [EC:", "\t[EC:", 1), "\t", 2) )
			} else { fds[11] = "" }

			TSV.Write ([]byte (strings.Join ( []string {
			fds[0], "-", "-", fds[3], fds[4], "-", "-", "-", "-", 
			fds[9], fds[10], fds[11]}, "\t") + "\n"))

			if isPathway {
				g := ""

				if strings.Contains (fds[4], "; ") {
					g = strings.Split (fds[4], "; ")[0]
				}

				G2P.Write ([]byte (fds[3] + "\t" + 
				strings.TrimLeft (fds[1], "PATH:") + "+" + fds[9] + "\t" + 
				strings.TrimLeft (strings.TrimRight (fds[11], "]"), "[EC:") + 
				"\t" + g + "\n") )
			}

		default:
			continue

		}
	}
}


func Match (name, datatsv string) (record []string, ok bool) {
	file, err := os.Open (datatsv)
	if err != nil { log.Fatal(err) }
	defer file.Close()

	species := formatSpeciesName (name)
	scanner := bufio.NewScanner (file)
	scanner.Scan () // skip header

	for scanner.Scan () {
		record = strings.Split (scanner.Text(), "\t")

		ok = (name == record[1] ||
		species == strings.Split (record[2], " (")[0])

		if ok { return }
	}

	return
}


func Update (saveto string) {
	log.Printf ("Update %s...\n", saveto)

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


func Get_local (code, path string) (ok bool) {
	Cmd := exec.Command ("cp", path + "/" + code + "00001.keg.gz", "./")
	err := Cmd.Run()

	if err != nil {
		log.Printf ("Failed to get %s from %s\n", code, path)
	} else { ok = true}

	return ok
}

func Get (codes []string) {
	ch := make (chan struct {}, 10)
	var wg sync.WaitGroup

	log.Printf ("Request organism code(s): %s\n", strings.Join (codes, " "))

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
	log.Printf ("Saved %s.gz\n", p)
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
