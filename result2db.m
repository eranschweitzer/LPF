function flag = result2db(tablename,colnames,T, varargin)

%% create connection
% since we are adding one row at a time, we leave the autocommit setting on
conn = varargin_parse(varargin,'conn', dbconn());

%% insert data
pause('on')
flag = true;
if ~isstruct(T)
	try
	    pause(randi(10)*1e-2)
	    datainsert(conn,tablename,colnames,vertcat(T{:}));
	catch
	    try
	        pause(randi(10)*1e-2)
	        datainsert(conn,tablename,colnames,vertcat(T{:}));
        catch 
	        flag = false;
	    end
	end
else
	for f = fieldnames(T).'
		flag_tmp = result2db(tablename.(f{1}), colnames.(f{1}), {T.(f{1})}, 'conn', conn);
		flag = flag && flag_tmp;
	end
end

%% close connection
if ~any(strcmp(varargin, 'conn'))
    close(conn)
end
