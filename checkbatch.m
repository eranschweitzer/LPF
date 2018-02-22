function flag = checkbatch(tablename,casename,IDbatch)

%% create connection
% since we are adding one row at a time, we leave the autocommit setting on
conn = dbconn();
setdbprefs('DataReturnFormat','numeric')
%% 
sql = ['SELECT count(*) FROM ', tablename, ' WHERE casename=''',casename,''''...
    ' and id=''', ids2str(unpackIDs(IDbatch(1,:))), ''''];
cnt = fetch(conn,sql);
flag = cnt ~= 0;
%% close connection
close(conn)