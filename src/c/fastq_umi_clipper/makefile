CC = gcc
CFLAGS = -O2
LIBDIR = /path/to/bstring/0.1.1/lib/
INCDIR = /path/to/bstring/0.1.1/include/bstring/
LIBS = -lbstring

	
fastq_umi_clipper:	fastq_umi_clipper.c 
	gcc -o fastq_umi_clipper -I${INCDIR} -L${LIBDIR} -lbstring fastq_umi_clipper.c 

clean:
	rm -f *.o fastq_umi_clipper 
