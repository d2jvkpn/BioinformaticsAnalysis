package main

import (
	"encoding/json"
	"encoding/xml"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"regexp"
	"strings"
	"sync"
)

const HELP = `
Query specices taxonomy ID, usage:
  $ biocrawler_TaxonID  Homo+sapiens  "Rattus norvegicus"  "Mus musculus"

author: d2jvkpn
version: 0.0.4
release: 2019-06-11
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main() {
	if len(os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
		fmt.Println(HELP)
		return
	}

	ch := make(chan struct{}, 10)
	var wg sync.WaitGroup

	for _, v := range os.Args[1:] {
		ch <- struct{}{}
		wg.Add(1)
		go func(p string, ch <-chan struct{}, wg *sync.WaitGroup) {
			defer func() { <-ch }()
			defer wg.Done()
			QueryTaxonId(p)
		}(formatSpeciesName(v), ch, &wg)
	}

	wg.Wait()

}

func QueryTaxonId(species string) {
	url := fmt.Sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"+
		"esearch.fcgi?db=taxonomy&term=\"%s\"",
		strings.Replace(species, " ", "+", -1))

	var err error
	var resp *http.Response
	var w Result
	var r *R1
	var bts []byte

	if resp, err = http.Get(url); err != nil {
		log.Printf("query failed of \"%s\": %v\n", species, err)
		return
	}

	defer resp.Body.Close()

	if bts, err = ioutil.ReadAll(resp.Body); err != nil {
		log.Printf("failed to read response for \"%s\": %v\n", species, err)
		return
	}

	if err = xml.Unmarshal(bts, &w); err != nil {
		log.Printf("xml unmarshal failed of \"%s\": %v\n", species, err)
		return
	}

	r = new(R1)
	r.Query = species
	r.ID = make([]int, 0, len(w.IdList))

	for i := range w.IdList {
		r.ID = append(r.ID, w.IdList[i].Id)
	}

	if bts, err = json.MarshalIndent(r, "", "    "); err != nil {
		log.Printf("json marshal failed of \"%s\": %v\n", species, err)
		return
	}

	bts = append(bts, []byte("\n\n")...)

	fmt.Println(string(bts))
}

type IdList struct{ Id int }
type Result struct{ IdList []IdList }

type R1 struct {
	Query string
	ID    []int
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
