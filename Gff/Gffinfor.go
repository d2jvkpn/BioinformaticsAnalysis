package main

import (
    "os"
    "fmt"
    "strings"
    "strconv"
    "bufio"
    "sort"
    "log"
    "text/tabwriter"
    "reflect"
    gzip "github.com/klauspost/pgzip"
)

const HELP = `
Summary gff/gtf (.gz) and extract attributions, usage: 
    summary sequences, sources, types
    $ Gffinfor  <gff>

    summary types' attributions
    $ Gffinfor  <gff>  <type1,type2...>

    extract attributions and "Dbxref" (tsv format)
    $ Gffinfor  <gff>  <type1,type2...>  <attr1,attr2,attr3...>  [dbxref1,dbxref2...]

author: d2jvkpn
version: 0.2
release: 2018-09-02
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisence: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

var parseAttr func (string, map[string]string)

func main () {
    if HasElem ([]string {"-h", "--help"}, os.Args[1]) { fmt.Println (HELP) }

    if strings.HasSuffix (os.Args[1], ".gtf") ||
        strings.HasSuffix (os.Args[1], ".gtf.gz") {
        parseAttr = gtfattr
    } else { parseAttr = gffattr }

    scanner, F, err := ReadInput (os.Args [1])
    if err != nil { log.Fatal (err) }
    if F != nil { defer F.Close () }

    switch len (os.Args) - 1 {
        case 1:
            P1 (scanner)
        case 2:
            P2 (scanner, strings.SplitN (os.Args[2], ",", -1) )
        case 3:
            P3 (scanner, strings.SplitN (os.Args[2], ",", -1), 
            strings.SplitN (os.Args[3], ",", -1) )
        case 4:
            P4 (scanner, strings.SplitN (os.Args[2], ",", -1), 
            strings.SplitN (os.Args[3], ",", -1),
            strings.SplitN (os.Args[4], ",", -1) )
        default:
            fmt.Println (HELP)
    }

}


//
func ReadInput (I string) (scanner *bufio.Scanner, F *os.File, err error) {
    if I == "-" {
        scanner = bufio.NewScanner (os.Stdin)
        return scanner, F, err
    }

    F, err = os.Open (I)
    if err != nil { return scanner, nil, err } 

    if strings.HasSuffix (I, ".gz") {
        gz, err := gzip.NewReader (F)
        if err != nil { return scanner, F, err }

        scanner = bufio.NewScanner (gz)
    } else {
        scanner = bufio.NewScanner (F)
    }

    return scanner, F, err
}


//
func TabPrint (array [][]string, head string) { 
    w := tabwriter.NewWriter (os.Stdout, 2, 0, 4, ' ', tabwriter.StripEscape)
    fmt.Fprintln (w, head + strings.Join (array[0], "\t"))
    for _, r := range array[1:] { fmt.Fprintln (w, head + strings.Join (r, "\t")) }
    w.Flush ()
}


func HasElem (s interface{}, elem interface{}) bool {
    arrV := reflect.ValueOf(s)

    if arrV.Kind() == reflect.Slice {
        for i := 0; i < arrV.Len(); i++ {
            // XXX - panics if slice element points to an unexported struct field
            // see https://golang.org/pkg/reflect/#Value.Interface
            if arrV.Index(i).Interface() == elem { return true }
        }
    }

    return false
}

//
func gtfattr (s string, kv map[string]string) {
    for _, i := range strings.Split (strings.TrimRight (s, ";"), "; ") {
        ii := strings.SplitN (i, " ", 2)
        kv [ii [0]] = strings.Trim (ii[1], "\"")
    }
}

func gffattr (s string, kv map[string]string) {
    for _, i := range strings.Split (s, ";") {
        ii := strings.SplitN (i, "=", 2)
        kv [ii [0]] = ii[1]
    }
}

//
func P1 (scanner *bufio.Scanner) {
    Sequences := make (map [string] int)
    Sources := make (map [string] int)
    Types := make (map [string] int)
    var fds []string

    for scanner.Scan() {
        line := scanner.Text ()
        if strings.HasPrefix (line, "#") { continue }
        fds = strings.SplitN (line, "\t", 9)
        Sequences[fds[0]] ++ 
        Sources[fds[1]] ++
        Types[fds[2]] ++
    }

    var array [][]string
    array = append (array, []string {"NAME", "COUNT"})

    array = append (array,
    []string {"Sequences", strconv.Itoa (len(Sequences))})

    var sKeys []string
    for k, _ := range Sources {	sKeys = append (sKeys, k) }
    // sort.Strings (sKeys)
    sort.Slice (sKeys, 
        func(i, j int) bool {
            return strings.ToLower (sKeys[i]) < strings.ToLower(sKeys[j]) 
        })

    for _, k:= range sKeys {
        x := []string {"source: " + k, strconv.Itoa (Sources[k])}
        array = append (array, x)
    }

    var tKeys []string
    for k, _ := range Types {tKeys = append (tKeys, k) }
    //sort.Strings (tKeys)
    sort.Slice (tKeys,
        func(i, j int) bool {
            return strings.ToLower (tKeys[i]) < strings.ToLower(tKeys[j])
        })

    for _, k:= range tKeys {
        x := []string {"type: " + k, strconv.Itoa (Types[k])}
        array = append (array, x)
    }

    TabPrint (array, "    ")
}

//
func P2 (scanner *bufio.Scanner, types []string) {

    TypeAttrs := make (map [string] map [string] int)
    var fds []string
    var array [][]string
    array = append (array, [] string {"TYPE\tATTRIBUTION", "TOTAL", "COUNT"})

    for scanner.Scan() {
        line := scanner.Text ()
        if strings.HasPrefix (line, "#") { continue }
        fds = strings.SplitN (line, "\t", 9)
        if types[0] != "" && ! HasElem (types, fds[2]) { continue }

        kv := make (map[string]string)
        parseAttr (fds[8], kv)

        for k, v := range kv {
           if _, ok := TypeAttrs [fds[2] + "\t" + k]; !ok {
                m := make(map[string] int)
                TypeAttrs [fds[2] + "\t" + k] = m
           }

           TypeAttrs [fds[2] + "\t" + k] [v] ++
        }
    }

    var keys []string
    for k, _ := range TypeAttrs { keys = append (keys, k) }

    sort.Slice (keys,
        func(i, j int) bool {
            return strings.ToLower (keys[i]) < strings.ToLower(keys[j])
        })

    for _, v := range keys {
        u := 0
        for _, c := range TypeAttrs[v] { u += c }

        array = append (array, 
        []string {v, strconv.Itoa (u), strconv.Itoa (len (TypeAttrs[v]))} )
    }

    TabPrint (array, "    ")
}

//
func P3 (scanner *bufio.Scanner, types []string, attrs []string) {
    var fds []string
    fmt.Println ( strings.Join (attrs, "\t") )

    for scanner.Scan() {
        line := scanner.Text ()
        if strings.HasPrefix (line, "#") { continue }
        fds = strings.SplitN (line, "\t", 9)
        if types[0] != "" && ! HasElem (types, fds[2]) { continue }

        kv := make (map[string]string)
        parseAttr (fds[8], kv)

        values := []string {}
        for _, k := range attrs { values = append (values, kv[k]) }

        fmt.Println (strings.Join (values, "\t"))
    }
}

//
func P4 (scanner *bufio.Scanner, types []string, attrs []string, dbx []string) {

    var fds []string
    fmt.Println (strings.Join (attrs, "\t") + "\t" + strings.Join (dbx, "\t") )

    for scanner.Scan() {
        line := scanner.Text ()
        if strings.HasPrefix (line, "#") { continue }
        fds = strings.SplitN (line, "\t", 9)
        if types[0] != "" && ! HasElem (types, fds[2]) { continue }

        kv := make (map[string]string)
        parseAttr (fds[8], kv)

        dkv := make (map[string]string)
        for _, d := range strings.Split (kv["Dbxref"], ",") {
            x := strings.SplitN (d, ":", 2)
            dkv[x[0]] = x[1]
        }

        values := []string {}
        for _, k := range attrs { values = append (values, kv[k]) }
        for _, k := range dbx { values = append (values, dkv[k]) }

        fmt.Println (strings.Join (values, "\t"))
    }
}
