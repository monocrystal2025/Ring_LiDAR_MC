function value = hash(input)
%HASH SHA-256 for deterministic cache input records.
arguments
    input
end
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(jsonencode(input),'UTF-8'));
value = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
