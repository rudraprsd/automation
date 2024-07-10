from pymatgen.io.vasp.outputs import Vasprun



def band_gap(vasprun_file):
    vaspout = Vasprun(vasprun_file)
    bandstr = vaspout.eigenvalue_band_properties()
    print(bandstr)
   

band_gap("../calculations/1951087_Si/bands/wsoc/band/vasprun.xml")
