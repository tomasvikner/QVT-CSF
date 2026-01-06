
function [METADATA, CLVALS, branchList, CBFSEG, Planes, StructCS, CSFSEG, CSFROI] = unpack_struct(data_struct)

    METADATA = data_struct.METADATA; 
    CLVALS = data_struct.CLVALS;
    branchList = data_struct.branchList;
    CBFSEG = data_struct.CBFSEG;
    Planes = data_struct.Planes;
    StructCS = data_struct.StructCS;
    CSFSEG = data_struct.CSFSEG;
    CSFROI = data_struct.CSFROI; % TEMP: prob not needed

end