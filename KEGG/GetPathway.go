package main

import (
  "os"
  "fmt"
  "log"
  "io/ioutil"
  "sync"
  "strings"
  "compress/gzip"
  "net/http"
)

const url = "http://www.kegg.jp/kegg-bin/download_htext?htext=%s&format=htext&filedir="

func main() {
  if len(os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help"{
    fmt.Println("Get KEGG pathway keg file by provide organism code(s), e.g. hsa mmu.")
    fmt.Println("\nproject: https://github.com/d2jvkpn/BioinformaticsAnalysis")
    return
  }

  var wg sync.WaitGroup

  for _, v := range(os.Args[1:]) {
    wg.Add(1)

    go func (p string) {
      defer wg.Done()

      resp, err := http.Get(fmt.Sprintf(url, p))
      if err != nil { log.Println (err); return }
      defer resp.Body.Close()
      // resp.Status, resp.StatusCode, resp.Proto, resp.ContentLength
      // resp.TransferEncoding, resp.Uncompressed

      body, err := ioutil.ReadAll(resp.Body)
      if err != nil { log.Println (err); return }
      lines := strings.Split (string(body), "\n")

      _b := strings.HasPrefix(lines[len(lines)-2], "#Last updated:")
      if ! _b { log.Printf ("Failed to scrap %s.\n", p); return }

      file, err := os.Create(p + ".gz")
      if err != nil { log.Println(err); return }
      defer file.Close()

      gw := gzip.NewWriter(file)
      gw.Write (body)
      gw.Close()

    } (v + "00001.keg")
  }

  wg.Wait()
}
