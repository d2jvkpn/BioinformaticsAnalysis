package main

import (
    "bufio"
    "fmt"
    "os"
    "log"
    "strings"
    "strconv"
    "runtime"
    gzip "github.com/klauspost/pgzip" //"compress/gzip"
)


var USAGE string = `
Usage: FastqCount  <input.fastq>  [phred]
    output (tsv): Total reads  Total bases  N bases  Q20  Q30  GC
    phred default: 33
    note: "pigz -dc *.fastq.gz | FastqCount -" is recommended for gzipped file(s).

author: d2jvkpn
version: 0.8
release: 2018-08-30
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main () {
    if len (os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
        fmt.Println (USAGE); return
    }

    var fname string = os.Args[1]
    var phred int = 33
    var err error

    if len (os.Args) == 3 {
        phred, err = strconv.Atoi (os.Args[2])
        if err != nil { log.Fatal (err) }
    }

    var scanner *bufio.Scanner

    if fname == "-" {
        scanner = bufio.NewScanner (os.Stdin)
    } else {
        file, err := os.Open (fname)
        if err != nil { log.Fatal (err) }

        if strings.HasSuffix (fname, ".gz") {
            gz, err := gzip.NewReader (file)
            if err != nil { log.Fatal (err) }
            scanner = bufio.NewScanner (gz)
        } else { scanner = bufio.NewScanner (file) }

        defer file.Close ()
    }

    var Rc, Bc, Q20, Q30, GC, Nc int = 0, 0, 0, 0, 0, 0
    ch :=  make (chan [2]string, 1000)
    runtime.GOMAXPROCS (2)

    go func() {
        for {
            var blk [2]string // !
            scanner.Scan(); scanner.Scan()
            blk[0] = scanner.Text ()
            scanner.Scan()
            if ! scanner.Scan() { break }
            blk[1] = scanner.Text ()
            ch <- blk
        }

        close (ch)
    } ()

    for k := range ch {
        Rc ++
        Bc += len (k[0])
        Nc += strings.Count (k[0], "N")
        GC += (strings.Count (k[0], "G") + strings.Count (k[0], "C"))

        for _, c := range k[1] {
            if int (c) - phred >= 20 { Q20 ++ } else { continue }
            if int (c) - phred >= 30 { Q30 ++ }
        }
    }

    fmt.Println ("Total reads\tTotal bases\tN bases\tQ20\tQ30\tGC")

    fmt.Printf ("%.2fM\t%.2fG\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n", 
        float64 (Rc) / float64 (1E+6),
        float64 (Bc) / float64 (1E+9),
        float64 (Nc*100) / float64 (Bc),
        float64 (Q20*100) / float64 (Bc),
        float64 (Q30*100) / float64 (Bc),
        float64 (GC*100) / float64 (Bc))

    fmt.Printf ("%d\t%d\t%d\t%d\t%d\t%d\n", Rc, Bc, Nc, Q20, Q30, GC)
}
