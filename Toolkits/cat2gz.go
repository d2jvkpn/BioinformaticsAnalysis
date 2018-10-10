package main

import (
	"fmt"
	"log"
	"os"
	"bufio"
	"flag"
	"strings"
	"path/filepath"
	gzip "github.com/klauspost/pgzip"
)

const USEAGE = `
Concatenate stdin, text files, gzipped files to one gzipped file.
Usage:
  $ cat2gz  <-o output.gz>  [-p cpu_max] [-level compress_level] \
  <input1.fastq input2.fastq.gz>

  note: when input is -, read standard input.`

const LISENSE = `
author: d2jvkpn
version: 0.5
release: 2018-10-09
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

func main() {
	outgz := flag.String("o", "", "output gzipped file")
	cpunum := flag.Int("p", 4, "max cpu number to use")
	level := flag.Int("level", 6, "compress level")
	flag.Parse()

	inputs := flag.Args()

	flag.Usage = func() {
		fmt.Println(USEAGE)
		flag.PrintDefaults()
		fmt.Println(LISENSE)
	}

	if len(inputs) == 0 || *outgz == "" { flag.Usage(); os.Exit(1) }

	type Writer interface {
		Write(p []byte) (n int, err error)
	}

	var wt Writer
	var err error

	err = os.MkdirAll(filepath.Dir(*outgz), 0755)
	if err != nil { log.Fatal(err) }

	var fwt *os.File

	fwt, err = os.Create(*outgz)
	if err != nil { log.Fatal(err) }
	defer fwt.Close()

	gznw, err := gzip.NewWriterLevel(fwt, *level)
	if err != nil { log.Fatal(err) }
	gznw.SetConcurrency(100000, *cpunum)

	wt = gznw
	defer gznw.Close()

	for _, s := range inputs {
		scanner, frd, err := ReadInput(s)
		if err != nil { log.Fatal(err) }
		defer frd.Close()

		for scanner.Scan() {
			_, err = wt.Write([]byte(scanner.Text() + "\n"))
			if err != nil { log.Fatal(err) }
		}
	}
}

func ReadInput(s string) (scanner *bufio.Scanner, file *os.File, err error) {
	if s == "-" {
		scanner = bufio.NewScanner(os.Stdin)
		return
	}

	file, err = os.Open(s)
	if err != nil { return }

	if strings.HasSuffix(s, ".gz") {
		var gz *gzip.Reader
		gz, err = gzip.NewReader(file)
		if err != nil { return }
		scanner = bufio.NewScanner(gz)
	} else {
		scanner = bufio.NewScanner(file)
	}

	return
}
