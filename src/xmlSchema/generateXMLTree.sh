rm schema.cpp
rm schema.h
xsdcxx cxx-tree --std c++11 schema.xsd
mv schema.cxx schema.cpp
mv schema.hxx schema.h
