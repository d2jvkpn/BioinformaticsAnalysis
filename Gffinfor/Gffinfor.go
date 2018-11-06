package main

import (
	"bufio"
	"fmt"
	"log"
	"net/url"
	"os"
	"errors"
	"strconv"
	"strings"
	"github.com/d2jvkpn/gopkgs/cmdplus"
)

const HELP = `
Summary gff/gtf (.gz) and extract attributions, usage:
    summary sequences, sources, types
    $ Gffinfor  <gff>
    note: stdin ("-") will be treated as gff format

    summary types' attributions
    $ Gffinfor  <gff>  <type1,type2...>

    extract attributions and Dbxref (tsv format)
    $ Gffinfor  <gff>  <type1,type2...>  <attr1,attr2,dbxref1,dbxref2...>

author: d2jvkpn
version: 0.8
release: 2018-10-27
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

var parseAttr func(string, map[string]string) error

func main() {
	if len(os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] ==  "--help" {
		fmt.Println (HELP)
		os.Exit (2)
	}

	if strings.HasSuffix(os.Args[1], ".gtf") ||
		strings.HasSuffix(os.Args[1], ".gtf.gz") {
		parseAttr = gtfattr
	} else {
		parseAttr = gffattr
	}

	scanner, frd, err := cmdplus.ReadCmdInput(os.Args[1])
	if err != nil { log.Fatal(err) }
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
func gtfattr(s string, kv map[string]string) (err error) {
	tmp := make (map[string]string)

	for _, i := range strings.Split(strings.TrimRight(s, ";"), "; ") {
		if i == "" { continue }
		ii := strings.SplitN(i, " ", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		if len(ii) != 2 {
			err = errors.New(fmt.Sprintf ("failed to split %s", i))
			return
		}

		tmp[ii[0]] = strings.Trim(ii[1], "\"")
	}

	for k, v := range tmp { kv[k] = v }

	return
}

func gffattr(s string, kv map[string]string) (err error) {
	tmp := make (map[string]string)

	for _, i := range strings.Split(s, ";") {
		if i == "" { continue }
		ii := strings.SplitN(i, "=", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		if len(ii) != 2 {
			err = errors.New (fmt.Sprintf ("failed to split %s", i))
			return
		}

		tmp[ii[0]] = ii[1]
	}

	for k, v := range tmp { kv[k] = v }

	return
}

//
func P1(scanner *bufio.Scanner) {
	Sequences := make(map[string]int)
	Sources := make(map[string]int)
	Types := make(map[string]int)
	var fds []string

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") { continue }
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" { continue }
		Sequences[fds[0]]++
		Sources[fds[1]]++
		Types[fds[2]]++
	}

	var array [][]string
	array = append(array, []string{"NAME", "COUNT"})

	array = append(array,
		[]string{"Sequences", strconv.Itoa(len(Sequences))})

	var sKeys []string
	for k, _ := range Sources { sKeys = append(sKeys, k) }

	cmdplus.SortStringSlice(sKeys)

	for _, k := range sKeys {
		x := []string{"source: " + k, strconv.Itoa(Sources[k])}
		array = append(array, x)
	}

	var tKeys []string
	for k, _ := range Types { tKeys = append(tKeys, k) }
	cmdplus.SortStringSlice(tKeys)

	for _, k := range tKeys {
		x := []string{"type: " + k, strconv.Itoa(Types[k])}
		array = append(array, x)
	}

	cmdplus.PrintStringArray(array)
}

//
func P2(scanner *bufio.Scanner, types []string) {
	TypeAttrs := make(map[string]map[string]int)
	var fds []string
	var err error
	var i int

	for scanner.Scan() {
		i ++
		line := scanner.Text()
		if strings.HasPrefix(line, "#") { continue }
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" { continue }
		if types[0] != "" && !cmdplus.HasElem(types, fds[2]) { continue }

		kv := make(map[string]string)
		err = parseAttr(fds[8], kv)

		if err != nil {
			log.Printf ("failed to parse attributions at line %d:" +
				"\n    %s\n    %s\n\n",  i, err, fds[8])

			continue
		}

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
	for k, _ := range TypeAttrs { keys = append(keys, k) }
	cmdplus.SortStringSlice(keys)

	for _, v := range keys {
		u := 0
		for _, c := range TypeAttrs[v] { u += c }

		array = append(array,
			[]string{v, strconv.Itoa(u), strconv.Itoa(len(TypeAttrs[v]))})
	}

	cmdplus.PrintStringArray(array)
}


//
func P3(scanner *bufio.Scanner, types, attrs []string) {
	var fds []string
	fmt.Println(strings.Join(attrs, "\t"))
	var err error
	var i int

	for scanner.Scan() {
		i ++
		line := scanner.Text()
		if strings.HasPrefix(line, "#") { continue }
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" { continue }
		if types[0] != "" && !cmdplus.HasElem(types, fds[2]) { continue }

		kv := make(map[string]string)
		err = parseAttr(fds[8], kv)

		if err != nil {
			log.Printf ("failed to parse attributions at line %d:" +
				"\n    %s\n    %s\n\n",  i, err, fds[8])

			continue
		}

		for _, d := range strings.Split(kv["Dbxref"], ",") {
			if d == "" { continue }
			x := strings.SplitN(d, ":", 2)
			kv[x[0]] = x[1]
		}

		values := []string{}
		for _, k := range attrs { values = append(values, kv[k]) }

		fmt.Println(strings.Join(values, "\t"))
	}
}
