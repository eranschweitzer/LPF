CREATE SCHEMA IF NOT EXISTS linpf;

CREATE TABLE IF NOT EXISTS linpf.res(
	casename VARCHAR(50),
	id VARCHAR(13),
	prop VARCHAR(4),
	criteria VARCHAR(3),
	value FLOAT
);

CREATE TABLE IF NOT EXISTS linpf.nogw(
	casename VARCHAR(50),
	id VARCHAR(13),
	prop VARCHAR(4),
	criteria VARCHAR(3),
	value FLOAT
);

CREATE TABLE IF NOT EXISTS linpf.actest(
	casename VARCHAR(50),
	id VARCHAR(13),
	dc BOOLEAN,
	linpf BOOLEAN
);

CREATE TABLE IF NOT EXISTS linpf.actestnogw(
	casename VARCHAR(50),
	id VARCHAR(13),
	dc BOOLEAN,
	linpf BOOLEAN
);
--CREATE INDEX on linpf.res(id);
--CREATE INDEX on linpf.res(casename,id);
--CREATE INDEX on linpf.res(casename,id,prop);
--CREATE INDEX on linpf.res(casename,id,criteria);
--CREATE INDEX on linpf.res(prop);
--CREATE INDEX on linpf.res(criteria);
--CREATE INDEX on linpf.res(prop,criteria);
--CREATE INDEX on linpf.res(id,prop,criteria);
--CREATE UNIQUE INDEX on linpf.res(casename,id,prop,criteria);
--
----Added later
--CREATE INDEX on linpf.res(casename);
--CREATE INDEX on linpf.res(value);
