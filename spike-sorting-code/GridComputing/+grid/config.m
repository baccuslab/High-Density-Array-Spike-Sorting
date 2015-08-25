function c = config()

c.submit_host = 'bs-submit01';
c.home = fullfile('/net/bs-filesvr02/export/group/hierlemann/');
%c.home = fullfile('/links/grid/scratch/hierlemann_sortings');

if ismac
    c.home = fullfile('/Users', getenv('USER'), 'tmp' )
end

end


