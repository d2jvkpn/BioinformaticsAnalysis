package main

import (
    "os"
    "fmt"
    "log"
    "io/ioutil"
    "strings"
    "time"
    "sync"
    "compress/gzip"
    "net/http"
)

const URL = "http://www.kegg.jp/kegg-bin/download_htext?htext=%s&format=htext&filedir="

const HELP = `
Download KEGG pathway keg file by provide organism code(s), e.g. hsa mmu.
$ Pathway_download hsa mmu ath

author: d2jvkpn
version: 0.5
release: 2018-09-02
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main () {
    if len (os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
        fmt.Println (HELP); return
    }

    num := 10
    if num > len (os.Args) - 1 { num = len (os.Args) - 1 }
    ch := make (chan struct{}, num)

    var wg sync.WaitGroup

    log.Printf ("%s Request organism code(s):\n    %s\n", 
    time.Now ().Format ("-0700"), strings.Join (os.Args[1:], " "))

    for _, v := range (os.Args[1:]) {
        ch <- struct{}{}
        wg.Add (1)
        go Querykeg (v + "00001.keg", ch, &wg)
    }

    wg.Wait ()
}

func Querykeg (p string, ch <- chan struct{}, wg *sync.WaitGroup) {
    defer func () { <- ch }()
    defer wg.Done ()

    log.Printf ("Querying %s...\n", p)

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
    return
}
