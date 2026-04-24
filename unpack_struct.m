
function [METADATA, CLVALS, branchList, CBFSEG, Planes, StructCS, CSFSEG, CSFROI, segment, flowCSF] = unpack_struct(DS)

    METADATA = DS.METADATA;
    if ~isfield(METADATA, 'iframes') || isempty(METADATA.iframes)
        METADATA.iframes = 1000; % standard high-res phase grid (interpCoupling)
    end
    CLVALS = DS.CLVALS;
    branchList = DS.branchList;
    CBFSEG = DS.CBFSEG;
    Planes = DS.Planes;
    if isfield(DS, 'StructCS')
        StructCS = DS.StructCS;
    else
        StructCS = struct();
    end
    if isfield(DS, 'CSFSEG')
        CSFSEG = DS.CSFSEG;
    else
        CSFSEG = struct();
    end
    if isfield(DS, 'CSFROI')
        CSFROI = DS.CSFROI; % TEMP: prob not needed
    else
        CSFROI = struct();
    end
    if isfield(DS, 'segment')
        segment = DS.segment;
    else
        segment = [];
    end
    if isempty(fieldnames(StructCS)) && ~isempty(fieldnames(CSFSEG))
        StructCS = CSFSEG; % single-source persisted payload compatibility
    end
    if isfield(DS, 'flowCSF')
        flowCSF = DS.flowCSF;
    elseif isfield(CLVALS, 'flowCSF')
        flowCSF = CLVALS.flowCSF;
    else
        flowCSF = struct();
    end

end