module gff3-go/gff3idx

go 1.25.0

require (
	gff3-go v0.0.0-00010101000000-000000000000
	github.com/zeebo/xxh3 v1.1.0
	golang.org/x/sys v0.46.0
)

require github.com/klauspost/cpuid/v2 v2.2.10 // indirect

replace gff3-go => ../
