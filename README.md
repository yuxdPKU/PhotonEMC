# PhotonEMC
Single photon analysis by using EMC

- Compile analysis module
```
cd module
mkdir build
cd build
../autogen.sh --prefix=$MYINSTALL
make -j 4
make install
```

- Make DST list of MDC2
```
cd macro/runList
./grablist.sh
```

- Submit jobs in condorJob
**Do not forget to change Initialdir**

- When job running, output root file temporarily stored in inReconstruction/[runnumber]

- When job finished, output root file saved in Reconstructed/[runnumber]
