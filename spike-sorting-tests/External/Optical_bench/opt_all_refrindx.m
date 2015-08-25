function n = opt_all_refrindx(wavelength)
% OPT_REFRINDX refractive index calculations.
% GLASS - name of glass type, WAVELENGTH in meters.
% Data from JML Optical Industries, Inc available at:
%   http://www.netacc.net/~jmlopt/transmission2.html
%   Or at:
%   http://www.jmlopt.com/level2/TechInfo/materialstable.aspx
% and
%
% Calling:
% n = opt_refrindx(glass,wavelength)
% 

% Copyright B. Gustavsson 20050804

persistent glass_names_sellmeier refr_consts_sellmeier

qwe = '';
fp = fopen('Sellmeier.glass.refr','r');
while ~feof(fp)
  qwe = str2mat(qwe,fgetl(fp));
end
glass_names_sellmeier = qwe(3:end,1:7);
refr_consts_sellmeier = str2num(qwe(3:end,8:end));


lambda_ref = [5876 4861 6563]*1e-10;

%I = strmatch(glass,glass_names_sellmeier);
%if ~isempty(I)
for I = 1:length(glass_names_sellmeier)  
  lambda = wavelength*1e6; % change units to micrometer
  A = refr_consts_sellmeier(I,1:3);
  B = refr_consts_sellmeier(I,4:6);
  n(I,:) = sqrt(1+A(1)*lambda.^2./(lambda.^2-B(1)) + A(2)*lambda.^2./(lambda.^2-B(2))+A(3)*lambda.^2./(lambda.^2-B(3)));
end

glasses = {'air','acrylic','b270','bak1','bak2','bak4','balkn3', ...
           'bk7','f2','f3','f4','fusedsilica','k5','k7','lasfn9', ...
           'lah71','pyrex','sapphire','sf1','sf2','sf8','sf10', ...
           'sf11','sf12','sf15','sf18','sf19','sf56','sk3','sk5', ...
           'sk11','sk16','ssk2','ssk4a','ssk51','zk5'};

for i = 1:length(glasses)
  glass = glasses{i};
  switch lower(glass)
   case 'air'
    nref = [1 1 1];
   case 'acrylic'
    nref = [1.491 1.496 1.488];
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
   otherwise
    lambda = [0 10000]*1e-9;
    alpha = [0 0];
    warning(['No values for glass absorption for: ',glass])
  end
  [abc,S,Mu] = polyfit(1./lambda_ref.^2,nref,2);
  n(I+i,:) = polyval(abc,1./wavelength.^2,[],Mu);
  
end

