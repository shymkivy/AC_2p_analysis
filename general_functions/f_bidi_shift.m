
function [data2, bidi_phase] = f_bidi_shift(data, n_frames, bidi_phase_frames)

[Ly, Lx, NT] = size(data);
data = reshape(data, Ly, Lx, 1, []);

if ~exist('n_frames', 'var'); n_frames = NT; end
if ~exist('n_shift_frames', 'var'); bidi_phase_frames = 1; end    

bidi_phase = zeros(n_frames,1);
data2 = zeros(size(data),'uint16');

disp('Computing offsets...');

if bidi_phase_frames > 1
    for yy = 1:n_frames-bidi_phase_frames
        bidi_phase(yy) = BiDiPhaseOffsets(data(:,:,1,yy:yy+bidi_phase_frames));
    end
    for yy = n_frames-bidi_phase_frames+1:n_frames
        bidi_phase(yy) = BiDiPhaseOffsets(data(:,:,1,yy:min(yy+bidi_phase_frames,NT)));
    end
else
    for yy = 1:n_frames-bidi_phase_frames
        bidi_phase(yy) = BiDiPhaseOffsets(data(:,:,1,yy));
    end
end

disp('Shifting...')
for yy = 1:n_frames
    data2(:,:,:,yy) = ShiftBiDi(bidi_phase(yy), data(:,:,:,yy), Ly, Lx);
end

data2 = reshape(data2, Ly, Lx, []);
disp('Done');
end