using Wavelets
using HDF5, JLD

type TestJLD
  Wvlt::Array
end
myvar = TestJLD([WT.db1])

JLD.save("test_file.jld","Wvlt", myvar.Wvlt)
tmp = JLD.load("test_file.jld")
