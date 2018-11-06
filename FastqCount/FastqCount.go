package main

import (
	"fmt"
	"log"
	"os"
	"flag"
	"path/filepath"
	"strings"
	"github.com/d2jvkpn/gopkgs/cmdplus"
)

const USAGE = `
Usage: FastqCount  [-phred value]  [-o tsv]  <input1.fastq input2.fastq.gz>
  output (tsv) header: Total reads  Total bases  N bases  Q20  Q30  GC
  note:
    1. When input is -, read standard input;
    2. "pigz -dc *.fastq.gz | FastqCount -" is recommended for gzipped file(s).
`

const LISENSE = `
author: d2jvkpn
version: 0.9.2
release: 2018-10-27
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main() {
	output := flag.String("o", "", "output summary to a tsv file, default: stdout")
	phred := flag.Int("phred", 33, "set phred value")

	flag.Usage = func() {
		fmt.Println (USAGE)
		flag.PrintDefaults ()
		fmt.Println (LISENSE)
	}

	flag.Parse()
	inputs := flag.Args()

	if len (os.Args) == 1 {
		flag.Usage()
		os.Exit(2)
	}

	type Writer interface { Write(p []byte) (n int, err error) }
	var err error
	var wt Writer
	ch := make(chan [2]string, 10000)
	var Rc, Bc, Q20, Q30, GC, Nc int = 0, 0, 0, 0, 0, 0

	go func() {
		for _, s := range inputs {
			log.Printf("FastqCount read sequeces from %s\n", s)

			scanner, frd, err := cmdplus.ReadCmdInput(s)
			if err != nil { log.Fatal(err) }
			defer frd.Close()

			for {
				var blk [2]string
				scanner.Scan()
				scanner.Scan()
				blk[0] = scanner.Text()
				scanner.Scan()
				if !scanner.Scan() { break }
				blk[1] = scanner.Text()
				ch <- blk
			}
		}
		close(ch)
	}()


	for k := range ch {
		Rc++
		Bc += len(k[0])
		Nc += strings.Count(k[0], "N")
		GC += strings.Count(k[0], "GC")

		for _, q := range k[1] {
			if int(q)-*phred >= 20 { Q20++ } else { continue }
			if int(q)-*phred >= 30 { Q30++ }
		}
	}

	wt = os.Stdout

	if *output != "" {
		err = os.MkdirAll (filepath.Dir (*output), 0755)
		if err != nil {
			log.Println (err)
		} else {
			wt, err = os.Create (*output)
			if err != nil { log.Println(err); wt = os.Stdout }
		}
	}

	fmt.Fprintln (wt, "Total reads\tTotal bases\tN bases\tQ20\tQ30\tGC")

	fmt.Fprintf (wt, "%.2fM\t%.2fG\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",
		float64(Rc) / float64(1E+6),
		float64(Bc) / float64(1E+9),
		float64(Nc*100) / float64(Bc),
		float64(Q20*100) / float64(Bc),
		float64(Q30*100) / float64(Bc),
		float64(GC*100) / float64(Bc))

	fmt.Fprintf (wt, "%d\t%d\t%d\t%d\t%d\t%d\n", Rc, Bc, Nc, Q20, Q30, GC)

	if *output != "" {
		log.Printf("Saved FastqCount result to %s\n", *output)
	}
}
