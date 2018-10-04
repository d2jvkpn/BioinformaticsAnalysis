package main

import (
	"bufio"
	"fmt"
	gzip "github.com/klauspost/pgzip"
	"log"
	"net/url"
	"os"
	"reflect"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
)

const HELP = `
Summary gff/gtf (.gz) and extract attributions, usage: 
    summary sequences, sources, types
    $ Gffinfor  <gff>

    summary types' attributions
    $ Gffinfor  <gff>  <type1,type2...>

    extract attributions and "Dbxref" (tsv format)
    $ Gffinfor  <gff>  <type1,type2...>  <attr1,attr2,dbxref1,dbxref2...>

author: d2jvkpn
version: 0.6
release: 2018-09-30
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

var parseAttr func(string, map[string]string)

func main() {
	if len(os.Args) == 1 || HasElem([]string{"-h", "--help"}, os.Args[1]) {
		fmt.Println(HELP)
		return
	}

	if strings.HasSuffix(os.Args[1], ".gtf") ||
		strings.HasSuffix(os.Args[1], ".gtf.gz") {
		parseAttr = gtfattr
	} else {
		parseAttr = gffattr
	}

	scanner, frd, err := ReadInput(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer frd.Close()

	switch len(os.Args) - 1 {
	case 1:
		P1(scanner)
	case 2:
		P2(scanner, strings.SplitN(os.Args[2], ",", -1))
	case 3:
		P3(scanner, strings.SplitN(os.Args[2], ",", -1),
			strings.SplitN(os.Args[3], ",", -1))
	default:
		fmt.Println(HELP)
		return
	}

}

//
func ReadInput(s string) (scanner *bufio.Scanner, file *os.File, err error) {
	if s == "-" {
		scanner = bufio.NewScanner(os.Stdin)
		return
	}

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

//
func TabPrint(array [][]string, leading string) {
	w := tabwriter.NewWriter(os.Stdout, 2, 0, 4, ' ', tabwriter.StripEscape)

	for _, r := range array {
		fmt.Fprintln(w, leading+strings.Join(r, "\t"))
	}

	w.Flush()
}

func HasElem(s interface{}, elem interface{}) bool {
	arrV := reflect.ValueOf(s)

	if arrV.Kind() == reflect.Slice {
		for i := 0; i < arrV.Len(); i++ {
			// XXX - panics if slice element points to an unexported struct field
			// see https://golang.org/pkg/reflect/#Value.Interface
			if arrV.Index(i).Interface() == elem {
				return true
			}
		}
	}

	return false
}

//
func gtfattr(s string, kv map[string]string) {
	for _, i := range strings.Split(strings.TrimRight(s, ";"), "; ") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, " ", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		kv[ii[0]] = strings.Trim(ii[1], "\"")
	}
}

func gffattr(s string, kv map[string]string) {
	for _, i := range strings.Split(s, ";") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, "=", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		kv[ii[0]] = ii[1]
	}
}

//
func P1(scanner *bufio.Scanner) {
	Sequences := make(map[string]int)
	Sources := make(map[string]int)
	Types := make(map[string]int)
	var fds []string

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		Sequences[fds[0]]++
		Sources[fds[1]]++
		Types[fds[2]]++
	}

	var array [][]string
	array = append(array, []string{"NAME", "COUNT"})

	array = append(array,
		[]string{"Sequences", strconv.Itoa(len(Sequences))})

	var sKeys []string
	for k, _ := range Sources {
		sKeys = append(sKeys, k)
	}

	SortStringSlice(sKeys)

	for _, k := range sKeys {
		x := []string{"source: " + k, strconv.Itoa(Sources[k])}
		array = append(array, x)
	}

	var tKeys []string
	for k, _ := range Types {
		tKeys = append(tKeys, k)
	}
	SortStringSlice(tKeys)

	for _, k := range tKeys {
		x := []string{"type: " + k, strconv.Itoa(Types[k])}
		array = append(array, x)
	}

	TabPrint(array, "	")
}

//
func P2(scanner *bufio.Scanner, types []string) {
	TypeAttrs := make(map[string]map[string]int)
	var fds []string

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		parseAttr(fds[8], kv)

		for k, v := range kv {
			nk := fds[2] + "\t" + k

			if _, ok := TypeAttrs[nk]; !ok {
				TypeAttrs[nk] = make(map[string]int)
			} else {
				TypeAttrs[nk][v]++
			}
		}
	}

	var array [][]string
	array = append(array, []string{"TYPE\tATTRIBUTION", "TOTAL", "UNIQUE"})

	var keys []string
	for k, _ := range TypeAttrs {
		keys = append(keys, k)
	}
	SortStringSlice(keys)

	for _, v := range keys {
		u := 0
		for _, c := range TypeAttrs[v] {
			u += c
		}

		array = append(array,
			[]string{v, strconv.Itoa(u), strconv.Itoa(len(TypeAttrs[v]))})
	}

	TabPrint(array, "	")
}


//
func P3(scanner *bufio.Scanner, types, attrs []string) {
	var fds []string
	fmt.Println(strings.Join(attrs, "\t"))

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		parseAttr(fds[8], kv)

		for _, d := range strings.Split(kv["Dbxref"], ",") {
			if d == "" {
				continue
			}
			x := strings.SplitN(d, ":", 2)
			kv[x[0]] = x[1]
		}

		values := []string{}
		for _, k := range attrs {
			values = append(values, kv[k])
		}

		fmt.Println(strings.Join(values, "\t"))
	}
}

//
func SortStringSlice(s []string) {
	sort.Slice(s, func(i, j int) bool {
		return strings.ToLower(s[i]) < strings.ToLower(s[j])
	})
}
