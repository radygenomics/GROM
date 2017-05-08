all: ./src/GROM.c

	tar xvzf ./src/samtools-1.3.1.tar.gz -C ./src
	$(MAKE) -C ./src/samtools-1.3.1
	$(MAKE) -C ./src/samtools-1.3.1/htslib-1.3.1

	gcc ./src/GROM.c -O2 -fopenmp -static -Wall -L./src/samtools-1.3.1 -L./src/samtools-1.3.1/htslib-1.3.1 -lbam -lhts -lm -lz -lpthread -I./src/samtools-1.3.1 -I./src/samtools-1.3.1/htslib-1.3.1 -o ./bin/GROM

clean:
	$(RM) ./bin/GROM
	$(MAKE) -C ./src/samtools-1.3.1 clean
	$(MAKE) -C ./src/samtools-1.3.1/htslib-1.3.1 clean

