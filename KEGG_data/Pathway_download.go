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

const HELP = `Archive KEGG pathway keg file by provide organism code(s), e.g. hsa mmu.
project: https://github.com/d2jvkpn/BioinformaticsAnalysis`

func main () {
    if len (os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
        fmt.Println (HELP); return
    }

    num := 10
    if num > len (os.Args) - 1 { num = len (os.Args) - 1 }
    ch := make (chan struct{}, num)

    var wg sync.WaitGroup

    log.Printf ("%s  request organism code(s): %s\n", 
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

    log.Printf ("%s  Querying %s...\n", time.Now ().Format ("-0700"), p)

    resp, err := http.Get (fmt.Sprintf (URL, p))
    if err != nil { log.Println (err); return }
    defer resp.Body.Close ()

    body, err := ioutil.ReadAll (resp.Body)
    if err != nil { log.Println (err); return }

    lines := strings.Split (string (body), "\n")
    if ! strings.HasPrefix (lines[len (lines)-2], "#Last updated:") {
        log.Printf ("%s  Failed to get %s\n", time.Now ().Format ("-0700"), p)
        return
    }

    file, err := os.Create (p + ".gz")
    if err != nil { log.Println (err); return }
    defer file.Close ()

    gw := gzip.NewWriter (file)
    gw.Write (body)
    gw.Close ()

    log.Printf ("%s  Saved %s...\n", time.Now ().Format ("-0700"), p)
    return
}
