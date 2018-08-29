package main

import (
    "bufio"
    "fmt"
    "os"
    "log"
    "strings"
    "strconv"
    gzip "github.com/klauspost/pgzip" // "compress/gzip"
)


var USAGE string = `Usage: FastqCount  <input.fastq>  [phred]
    Phred default: 33
    Output (tsv): Total reads  Total bases  N bases  Q20  Q30  GC
    Note: "pigz -dc *.fastq.gz | FastqCount -" is recommended for gzipped file(s).

author: d2jvkpn
version: 0.7
release: 2018-08-28
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
        if err != nil { log.Fatal (err)}
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

        if err != nil { log.Fatal (err) }
        defer file.Close ()
    }

    var Rc, Bc, Q20, Q30, GC, Nc int = 0, 0, 0, 0, 0, 0
    var sequence, qual string

    for {
        scanner.Scan(); scanner.Scan()
        sequence = scanner.Text ()
        scanner.Scan()
        if ! scanner.Scan() { break }
        qual = scanner.Text ()

        Rc ++
        Bc += len(sequence)
        Nc += strings.Count (sequence, "N")
        GC += (strings.Count (sequence, "G") + strings.Count (sequence, "C"))

        for _, c := range qual {
            if int (c) - phred >= 20 { Q20 ++ } else { continue }
            if int (c) - phred >= 30 { Q30 ++ }
        }
    }

    fmt.Println ("Total reads\tTotal bases\tN bases\tQ20\tQ30\tGC")

    fmt.Printf ("%d (%.2fM)\t%d (%.2fG)\t%d (%.2f%%)\t%.2f%%\t%.2f%%\t%.2f%%\n", 
        Rc, float32 (Rc) / 1E+6,
        Bc, float32 (Bc) / 1E+9,
        Nc, float32 (Nc*100 / Bc),
        float32 (Q20*100 / Bc),
        float32 (Q30*100 / Bc),
        float32 (GC*100 / Bc))
}
