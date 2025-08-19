% Area adaptive local thresholding of the CSF segment 
function [bics, bcsf, cnob, bcube, full, stdthresh, thlist] = segment_pcsf(Slice, BICSTHRESH, STDTHRESH, segment)

    nvcbf = sum(segment(:));

    SE = strel('sphere', 20); % 8 voxels/mm 

    % make sure to exclude out of FoV 
    full = Slice.mcsf > 0; 
    slicemask = zeros(size(full));
    slicemask(2:end-1, 2:end-1) = 1; 
    full = full .* slicemask;

    b1f = 2; 
    b3f = 4; 
    thlist.bt1 = [0.5, 0.75, 1.0, 1.25, 1.5]; 
    thlist.bt3 = [0.5, 1, 1.5, 2, 2.5]; 
    nthresh = numel(thlist.bt1);
    thlist.bics = cell(nthresh, nthresh);
    thlist.bcsf = cell(nthresh, nthresh);
    thlist.cnob = cell(nthresh, nthresh);

    for b1 = 1:nthresh

        for b3 = 1:nthresh

            bt1 = thlist.bt1(b1);
            bt3 = thlist.bt3(b3);
    
            count = 0;
            nvcbics = 0; % num voxels after BICS thresholding 
            bicsthresh = BICSTHRESH;
            while (nvcbics < (bt1 * nvcbf)) && (count < 50) % if blood seg is more than twice as large, adapt threshold 
                bics = (Slice.mcsf < bicsthresh) & full;
                cc = bwconncomp(bics);
                nc = cc.NumObjects;
                if nc > 1
                    [rows, cols] = size(bics);
                    center = [rows, cols] / 2;
                    stats = regionprops(cc, 'Centroid');
                    dists = zeros(nc, 1);
                    for i = 1:nc
                        c = stats(i).Centroid;
                        dists(i) = norm(c - center);
                    end
                    [~, idx] = min(dists);
                    central_bics = false(size(bics));
                    central_bics(cc.PixelIdxList{idx}) = true;
                    cbics = central_bics;
                else
                    cbics = bics;
                end
                nvcbics = sum(cbics(:)); % num voxels in central bics 
                bicsthresh = bicsthresh + 0.01; 
                count = count + 1; 
            end
            
            % If central BICS is too large 
            count = 0; 
            while (nvcbics > (bt1 * nvcbf)) && (count < 50) % if blood seg is more than twice as large, adapt threshold 
                bics = (Slice.mcsf < bicsthresh) & full;
                cc = bwconncomp(bics);
                nc = cc.NumObjects;
                if nc > 1
                    [rows, cols] = size(bics);
                    center = [rows, cols] / 2;
                    stats = regionprops(cc, 'Centroid');
                    dists = zeros(nc, 1);
                    for i = 1:nc
                        c = stats(i).Centroid;
                        dists(i) = norm(c - center);
                    end
                    [~, idx] = min(dists);
                    central_bics = false(size(bics));
                    central_bics(cc.PixelIdxList{idx}) = true;
                    cbics = central_bics;
                else
                    cbics = bics;
                end
                nvcbics = sum(cbics(:)); % num voxels in central bics 
                bicsthresh = bicsthresh - 0.01; 
                count = count + 1; 
            end
        
            dcbics = imdilate(cbics, SE); 
        
            count = 0;
            nvstd = 0; % num voxels after SD thresholding 
            temp = Slice.scsf .* dcbics;
            stdthresh = STDTHRESH;
            while (nvstd < (bt3 * nvcbf)) && (count < 100) % should always enter at least once 
                sdflow = temp > stdthresh;
                pacsf = sdflow & ~bics;
                nvstd = sum(pacsf(:));
                stdthresh = stdthresh - 0.01; 
                count = count + 1; 
            end
            % adapt threshold if CSF segment is too large
            count = 0; 
            while (nvstd > (bt3 * nvcbf)) && (count < 100)
                sdflow = temp > stdthresh;
                pacsf = sdflow & ~bics;
                nvstd = sum(pacsf(:));
                stdthresh = stdthresh + 0.01; 
                count = count + 1; 
            end
            sdflow = temp > stdthresh;
            pacsf = sdflow & ~bics;
        
            if nvstd > (bt3 * nvcbf)
                disp(' *** COUNT TO ALTER CSF SD THRESH INSUFFUCIENT. *** ')
            end
        
            bics = bics & full;
            bcsf = sdflow & full; 
            cnob = pacsf & full; 
            bcube = Slice.cube > graythresh(Slice.cube) & full; 

            thlist.bics{b1, b3} = bics;
            thlist.bcsf{b1, b3} = bcsf;
            thlist.cnob{b1, b3} = cnob;

        end

    end

    % make sure to return to standard ones 
    bics = thlist.bics{b1f, b3f};
    bcsf = thlist.bcsf{b1f, b3f};
    cnob = thlist.cnob{b1f, b3f};

end