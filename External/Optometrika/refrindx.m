function n_refr = refrindx( wavelength, glass )
% REFRINDX refractive index calculations.
% Data from JML Optical Industries, Inc available at:
%   http://www.netacc.net/~jmlopt/transmission2.html
%   Or at:
%   http://www.jmlopt.com/level2/TechInfo/materialstable.aspx
%
%   The human eye data are from in Escudero-Sanz & Navarro, 
%   "Off-axis aberrations of a wide-angle schematic eye model", 
%   JOSA A, 16(8), 1881-1891 (1999).

persistent glass_names_sellmeier refr_consts_sellmeier

if isempty( glass_names_sellmeier )
    qwe = '';
    fp = fopen( 'Sellmeier.glass.refr', 'r' );
    while ~feof( fp )
        qwe = str2mat( qwe, fgetl( fp ) );
    end
    glass_names_sellmeier = qwe( 2:end, 1:7 );
    refr_consts_sellmeier = str2num( qwe( 2:end, 8:end ) );
end

lambda_ref = [ 5876 4861 6563 ] * 1e-10;
lambda_eye = [4580 5430 5893 6328] * 1e-10;

ind = strmatch( glass, glass_names_sellmeier );
if ~isempty( ind )
    lambda = wavelength * 1e6; % change units to micrometer
    A = refr_consts_sellmeier( ind, 1:3 );
    B = refr_consts_sellmeier( ind, 4:6 );
    n_refr = sqrt( 1 + A(1) * lambda.^2 ./ ( lambda.^2 - B(1) ) + ...
        A(2) * lambda.^2 ./ ( lambda.^2 - B(2) ) + ...
        A(3) * lambda.^2 ./ ( lambda.^2 - B(3) ) );
else
    switch lower(glass)
        case 'air'
            nref = [1 1 1];
        % plastics
        case { 'pmma', 'acrylic' }
            nref = [1.491 1.496 1.488];
        case { 'pc', 'polycarbonate' }
            nref = [1.5849 1.5994 1.5782];
        case { 'ps', 'polystyrene' }
            nref = [1.5917 1.6056 1.5853];
        case { 'nas-21' 'nas' }
            nref = [1.5714 1.5835 1.5669];
        case { 'optorez-1330' 'optorez' }
            nref = [1.5094 1.5163 1.5067];
        case { 'zeonex-e48r' 'zeonex' }
            nref = [1.5309 1.5376 1.5273];
        % glasses
        case 'b270'
            nref = [1.5230 1.5292 1.5202];
        case 'bak1'
            nref = [1.5725 1.5794 1.5694];
        case 'bak2'
            nref = [1.5399 1.5462 1.5372];
        case 'bak4'
            nref = [1.56883 1.5759 1.56576];
        case 'balkn3'
            nref = [1.51849 1.52447 1.51586];
        case 'bk7'
            nref = [1.5168 1.5224 1.5143];
        case 'f2'
            nref = [1.6200 1.6320 1.6150];
        case 'f3'
            nref = [1.61293 1.62461 1.60806];
        case 'f4'
            nref = [1.6165 1.6284 1.6116];
        case 'fusedsilica'
            nref = [1.458 1.463 1.456];
        case 'k5'
            nref = [1.5224 1.5285 1.5198];
        case 'k7'
            nref = [1.51112 1.517 1.50854];
        case 'lasfn9'
            nref = [1.850 1.8689 1.8425];
        case 'lah71'
            nref = [1.8502 1.8689 1.8425];
        case 'pyrex'
            nref = [1.473 1.478 1.471];
        case 'sapphire'
            nref = [1.7682 1.7756 1.7649];
        case 'sf1'
            nref = [1.71736 1.73462 1.71031];
        case 'sf2'
            nref = [1.6476 1.6612 1.6421];
        case 'sf5'
            nref = [1.6727 1.6875 1.66661];
        case 'sf8'
            nref = [1.6889 1.7046 1.6825];
        case 'sf10'
            nref = [1.72825 1.74648 1.72085];
        case 'sf11'
            nref = [1.7847 1.8064 1.7759];
        case 'sf12'
            nref = [1.64831 1.66187 1.64271];
        case 'sf15'
            nref = [1.69895 1.71546 1.69221];
        case 'sf18'
            nref = [1.7215 1.7390 1.7143];
        case 'sf19'
            nref = [1.6668 1.6811 1.6609];
        case 'sf56'
            nref = [1.7847 1.8061 1.7760];
        case 'sk3'
            nref = [1.6088 1.6160 1.6056];
        case 'sk5'
            nref = [1.5891 1.5958 1.5861];
        case 'sk11'
            nref = [1.5638 1.5702 1.5610];
        case 'sk16'
            nref = [1.6204 1.6275 1.6172];
        case 'ssk2'
            nref = [1.6223 1.63048 1.61878];
        case 'ssk4a'
            nref = [1.61765 1.62547 1.61427];
        case 'ssk51'
            nref = [1.60361 1.61147 1.60022];
        case 'zk5'
            nref = [1.53375 1.54049 1.53084];
        % eye
        case 'cornea'
            nref = [1.3828 1.3777 1.376 1.3747];
        case 'aqueous'
            nref = [1.3445 1.3391 1.3374 1.336];
        case 'lens'
            nref = [1.4292 1.4222 1.42 1.4183];
        case 'vitreous'
            nref = [1.3428 1.3377 1.336 1.3347];
        otherwise
            error( [ 'No values for glass absorption for: ', glass ] );
    end
    if strcmp( glass, 'cornea' ) || ...
            strcmp( glass, 'aqueous' ) || ...
            strcmp( glass, 'lens' ) || ...
            strcmp( glass, 'vitreous' )
        [ abc, ~, Mu ] = polyfit( 1./lambda_eye.^2, nref, 2 );
    else
        [ abc, ~, Mu ] = polyfit( 1./lambda_ref.^2, nref, 2 );
    end
    n_refr = polyval( abc, 1./wavelength.^2, [], Mu );
end

end
