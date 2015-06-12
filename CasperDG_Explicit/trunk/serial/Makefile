
# HELP: displays information message
# ---------------------------------------------------
help    :
	@echo ""
	@echo " Usage:"
	@echo "	make [TARGET]"
	@echo ""
	@echo " TARGET:"
	@echo "	pgi32 ... builds binaries using pgi  "
	@echo "	pgi64 ... builds binaries using pgi64"
	@echo "	ifort ... builds binaries using ifort"
	@echo "	gfortran ... builds binaries using gfort"
	@echo "	g95   ... builds binaries using g95  "
	@echo""
	@echo "	help  ... displays this message      "
	@echo "	clean ... removes object files       "
	@echo "	backup ... creates tarball file       "

# Compiler Targets
# ----------------------------------------------
% :
	@echo "                                      "
	@echo " *** Compiling : $@ version           "
	@echo "                                      "
	$(MAKE) -f Makefile.dg ARCH=$@

# Maintanace Targets
# ----------------------------------------------
clean   :
	@echo "                                      "
	@echo " *** Clean : removing object files    "
	@echo "                                      "
	@( $(RM) *.o *.mod *~ )
	@( $(RM) ../src/*~ )

backup  :
	@echo "                                      "
	@echo " *** Backup : tarball files           "
	@echo "                                      "
	@( tar -zcf tarball.tar.gz *.F90 Makefile*)
