module github.com/EndCredits/gff3-go/gff3idx

go 1.22

require (
	github.com/EndCredits/gff3-go v0.0.0-00010101000000-000000000000
	github.com/zeebo/xxh3 v1.0.2
	golang.org/x/sys v0.22.0
)

require github.com/klauspost/cpuid/v2 v2.2.9 // indirect

replace github.com/EndCredits/gff3-go => ../
