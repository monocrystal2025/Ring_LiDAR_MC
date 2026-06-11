% Make manually launched MATLAB sessions discoverable by the MCP server.
try
    shareMATLABSession();
    fprintf('MATLAB MCP session sharing enabled.\n');
catch exception
    warning('startup:MATLABMCP', ...
        'Could not enable MATLAB MCP session sharing: %s', exception.message);
end
