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
KEGG pathway process, usage:

1. update local data table for command "match":
    $ Pathway  Update

2. find match species name or code in local data table:
    $ Pathway  match  "Rhinopithecus roxellana"
    $ Pathway  match  Rhinopithecus+roxellana
    $ Pathway  match  rro

3. get organisms keg file(s) by download:
    $ Pathway  Get  hsa mmu ath

4. download pathway html:
    $ Pathway  HTML  hsa00001.keg.gz  ./hsa00001
    Note: existing html files will not be overwritten

5. get organisms a single keg file from local:
    $ Pathway  get  hsa

6. convert keg format to tsv:
    $ Pathway  totsv  hsa00001.keg.gz  hsa00001.keg.tsv

    hsa00001.keg.tsv tsv header:
    C_id C_entry C_name gene_id gene A_id A_name B_id B_name KO KO_name EC

    other output files: gene_id2pathway.tsv pathway.infor.tsv

6. get species keg (from local), convert to tsv and download html files:
    $ Pathway  species  Rhinopithecus+roxellana
    output files: rro00001.keg.gz rro00001.keg.tsv

author: d2jvkpn
version: 0.8
release: 2018-09-14
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main () {
	nargs := len (os.Args) - 1

	if nargs == 0 || os.Args[1] == "-h" ||  os.Args[1] == "--help" {
		fmt.Println (HELP); return
	}

	cmd := os.Args[1]
	ep, _ := exec.LookPath (os.Args[0])
	datatsv := filepath.Dir (ep) + "/KEGG_data/KEGG_organism.tsv"

	switch {
	case nargs == 1 && cmd == "Update":
		Update (datatsv)

	case nargs > 1 && cmd == "Get":
		Get (os.Args[2:])

	case nargs == 2 && cmd == "get":
		ok := Get_local (os.Args[2], filepath.Dir (ep) + "/KEGG_data/Pathway_kegs")
		if !ok { os.Exit(1) }

	case nargs == 3 && cmd == "HTML":
		DownloadHTML (os.Args[2], os.Args[3])

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

			DownloadHTML (record[1] + "00001.keg.gz", record[1] + "00001")
		} else {
			fmt.Println ("NotFound")
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


func DownloadHTML (keg, outdir string) {
    var bts []byte
    var err error

    frd, err := os.Open (keg)
	if err != nil { log.Fatal (err) }
	defer frd.Close()

	err = os.MkdirAll (outdir, 0755)
	if err != nil { log.Fatal (err)}

	log.Printf ("Save html files to %s\n", outdir)

	if strings.HasSuffix (keg, ".gz") {
		var gz *gzip.Reader
		gz, err = gzip.NewReader (frd)
		if err != nil { return }

		bts, err = ioutil.ReadAll (gz)
	} else {
		bts, err = ioutil.ReadAll (frd)
	}

    re := regexp.MustCompile ("PATH:[a-z]+[0-9]+")
    PATHs := re.FindAllString (string (bts), -1)

	for i, _ := range (PATHs) { PATHs[i] = strings.TrimLeft (PATHs[i], "PATH:") }

    var wg sync.WaitGroup
	ch := make (chan struct {}, 10)

	for _, p := range PATHs {
		ch <- struct{}{}
		wg.Add (1)
		go gethtml (p, outdir, ch, &wg)
	}
    wg.Wait()
}


func gethtml (p, outdir string, ch <- chan struct{}, wg *sync.WaitGroup) {
	defer func () { <- ch }()
	defer wg.Done ()

	html := outdir + "/" + p + ".html"
	png := outdir + "/" + p + ".png"
	url := "http://www.genome.jp/kegg"
	code := p[:(len(p)-5)]

	if _, err := os.Stat (html); err == nil { return }

	htmlurl := fmt.Sprintf (url + "-bin/show_pathway?%s", p)
	pngurl := fmt.Sprintf (url + "/pathway/%s/%s.png", code, p)

	htmlresp, err := http.Get (htmlurl)
	if err != nil { log.Println (err); return }
	defer htmlresp.Body.Close ()

	htmlbody, err := ioutil.ReadAll (htmlresp.Body)
	if err != nil { log.Println (err); return }

	text := string (htmlbody)
	if ! strings.HasSuffix (text, "</html>\n") { return }

	re, _ := regexp.Compile("\\[[\\S\\s]+?\\]")  
	text = re.ReplaceAllString (text, "")

	re, _ = regexp.Compile("\\<script[\\S\\s]+?\\</script\\>")  
	text = re.ReplaceAllString (text, "")

	re, _ = regexp.Compile("\\<style[\\S\\s]+?\\</style\\>")  
	text = re.ReplaceAllString (text, "")

	re, _ = regexp.Compile("\\<link[\\S\\s]+?\\</link\\>")  
	text = re.ReplaceAllString (text, "")

	re, _ = regexp.Compile("\\<table[\\S\\s]+?\\</table\\>")  
	text = re.ReplaceAllString (text, "")

	text = strings.Replace (text, fmt.Sprintf("/kegg/pathway/%s/", code), "", 1)

	text = strings.Replace (text, "/dbget-bin/www_bget?", 
	"https://www.genome.jp/dbget-bin/www_bget?", -1)

	pngresp, err := http.Get (pngurl)
	if err != nil { log.Println (err); return }
	defer pngresp.Body.Close ()

	pngbody, err := ioutil.ReadAll (pngresp.Body)
	if err != nil { log.Println (err); return }

	err = ioutil.WriteFile (png, pngbody, 0664)
	if err != nil {log.Println (err); return}

	err = ioutil.WriteFile (html, []byte (text), 0664)
	if err != nil {log.Println (err); return}
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
	var fds [11]string
	var isPathway bool

	TSV.Write ([] byte ("C_id\tC_entry\tC_name\tgene_id\tgene\t" + 
	"A_id\tA_name\tB_id\tB_name\tKO\tEC\n"))

	G2P.Write ([]byte ("gene_id\tpathway\tKO\tEC\tgene_name\n"))
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

			fds[3], fds[4], fds[9], fds[10]= "", "", "", ""
			TSV.Write ([]byte (strings.Join (fds[0:], "\t")  + "\n"))
			isPathway = strings.HasPrefix (fds[1], "PATH:")

			if isPathway {
				Pinfor.Write ([]byte ( strings.TrimLeft (fds[1], "PATH:") +
				"\t" + fds[2] + "\t" + fds[6] + "\t" + fds[8] + "\n"))
			}

		case 'D':
			tmp := strings.SplitN (strings.TrimLeft (line, "D	  "), "\t", 2)
			if len(tmp) != 2 { continue }

			copy (fds[3:5], strings.SplitN (tmp[0], " ", 2))

			if strings.Contains (tmp[1], " [EC:") {
				copy (fds[9:11], strings.SplitN (
				strings.Replace (tmp[1], " [EC:", "\t[EC:", 1), "\t", 2) )
			} else { fds[10] = "" }

			TSV.Write ([]byte (strings.Join ( []string {
			fds[0], "-", "-", fds[3], fds[4], "-", "-", "-", "-", 
			fds[9], fds[10]}, "\t") + "\n"))

			if isPathway {
				g := ""
				KO := strings.Split (fds[9], " ")[0]

				if strings.Contains (fds[4], "; ") {
					g = strings.Split (fds[4], "; ")[0]
				}

				G2P.Write ([]byte (fds[3] + "\t" + 
				strings.TrimLeft (fds[1], "PATH:") + "\t" + KO + "\t" + 
				strings.TrimLeft (strings.TrimRight (fds[10], "]"), "[EC:") + 
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
	log.Printf ("Update %s\n", saveto)

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

	resp, err := http.Get (fmt.Sprintf ("http://www.kegg.jp/kegg-bin" + 
	"/download_htext?htext=%s&format=htext&filedir=", p))

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
