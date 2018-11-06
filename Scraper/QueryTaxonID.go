package main

import (
	"encoding/xml"
	"fmt"
	"io/ioutil"
	"net/http"
	"os"
	"regexp"
	"strings"
	"sync"
)

const HELP = `
Query specices taxonomy ID, usage:
    QueryTaxonID <species1> <species2> ...

    $ QueryTaxonId Homo+sapiens "Rattus norvegicus" "Mus musculus"
    Outputs:
    Rattus norvegicus: 10116
    Homo sapiens: 9606
    Mus musculus: 10090

author: d2jvkpn
version: 0.0.3
release: 2018-09-16
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
			defer wg.Done()
			defer func() { <-ch }()

			QueryTaxonId(p)
		}(formatSpeciesName(v), ch, &wg)
	}

	wg.Wait()

}

func QueryTaxonId(species string) {
	url := fmt.Sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"+
		"esearch.fcgi?db=taxonomy&term=\"%s\"",
		strings.Replace(species, " ", "+", -1))

	resp, err := http.Get(url)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: query failed, %s\n", species, err.Error())

		return
	}

	defer resp.Body.Close()

	body, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: query failed, %s\n", species, err.Error())

		return
	}

	var w Result
	err = xml.Unmarshal(body, &w)

	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: query failed, %s\n", species, err.Error())

		return
	} else {
		fmt.Printf(species + ": " + w.IdList[0].Id + "\n")
	}
}

type IdList struct{ Id string }

type Result struct{ IdList []IdList }

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
