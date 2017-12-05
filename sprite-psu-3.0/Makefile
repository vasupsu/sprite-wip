include Make.inc
all:
	@echo "Building PRUNE"
	cd prune; $(MAKE); cd ..;
#
	@echo "Building SAMPA and PARSNIP"
	cd parsnip; $(MAKE); cd ..;
#
	@echo "Building utility to view AEB,AIB files"
	cd utils; $(MAKE); cd ..;
#
	@echo "DONE"
clean:
	rm -f prune/*.o prune/prune-bwamem prune/prune-genreadidx prune/bam2aebaib \
	parsnip/parsnip parsnip/sampa utils/aebaib2txt
