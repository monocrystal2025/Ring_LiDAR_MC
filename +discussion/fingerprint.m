function digest = fingerprint()
%FINGERPRINT SHA-256 of the shared simulator source for cache validation.
folder = fileparts(mfilename('fullpath'));
files = dir(fullfile(folder,'*.m'));
[~,order] = sort(string({files.name}));
md = java.security.MessageDigest.getInstance('SHA-256');
for k = order
    contents = unicode2native([files(k).name,newline, ...
        fileread(fullfile(folder,files(k).name))],'UTF-8');
    md.update(contents);
end
bytes = typecast(md.digest(),'uint8');
digest = string(lower(reshape(dec2hex(bytes,2).',1,[])));
end
