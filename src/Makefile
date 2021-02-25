RM = rm
ifeq ($(OS),Windows_NT)
    RM = cmd.exe /C del
endif

FC = gfortran
BULK = bulkm.f90 bulkio.f90 bulk.f90 bulksub.f90
SURF = surfm.f90 surf.f90 surfio.f90 surfsub.f90 domainsum.f90
DOMAIN = domain.f90 domainsum.f90
XYZ = xyzb.f90 strfac.f90 asf.f90 
SOURCES = asf.f90 strfac.f90 scpot.f90 trmatg.f90 gcmi.f90
LAPACK = zgeev.f zgeevb.f

ALL: bulk surf

bulk: $(BULK) $(SOURCES)
	$(FC) -O3 -o bulk.exe $(BULK) $(SOURCES) $(LAPACK)

surf: $(SURF) $(SOURCES)
	$(FC) -O3 -o surf.exe $(SURF) $(SOURCES) $(LAPACK)

domain: $(DOMAIN)
	$(FC) -o domain.exe $(DOMAIN)

xyz: $(XYZ)
	$(FC) -o xyz.exe $(XYZ)

clean:
	-$(RM) bulk.exe surf.exe domain.exe xyz.exe

