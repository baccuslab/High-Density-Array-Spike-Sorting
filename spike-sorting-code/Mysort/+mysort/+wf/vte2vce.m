function vce = vte2vce(vte, nC)
    % inverse of vce2vte, see explanation there
    
    if isempty(vte)
        vce = [];
        return
    end
    
    Tf = size(vte,2)/nC;
    assert(Tf==round(Tf), 'nC does not match vce !');
    nT = size(vte,1);
    
    vce = reshape(permute(reshape(vte', [Tf nC nT]), [2 1 3]), [nC*Tf nT])';