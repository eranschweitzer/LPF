CREATE TABLE IF NOT EXISTS linpf.vt(
	caseid SMALLINT,
	id VARCHAR(12),
	prop VARCHAR(3),
	criteria VARCHAR(3),
	value FLOAT
);

COMMENT ON COLUMN linpf.vt.criteria IS 'All voltages in per unit, theta in radians, del in percent';

CREATE TABLE IF NOT EXISTS linpf.flows(
	caseid SMALLINT,
	pfid VARCHAR(12),
	id VARCHAR(5),
	prop VARCHAR(2),
	criteria VARCHAR(3),
	value FLOAT
);

COMMENT ON COLUMN linpf.flows.criteria IS 'All flows in p.u., del in percent';

CREATE TABLE IF NOT EXISTS linpf.cases(
	id SMALLSERIAL,
	name VARCHAR(256),
	PRIMARY KEY(id)
);

INSERT INTO linpf.cases (name) VALUES
	('case24_ieee_rts'),
	('case30'),
	('case57'),
	('case118'),
	('case145'),
	('case300'),
	('case1354pegase'),
	('case1888rte'),
	('case1951rte'),
	('case2383wp'),
	('case2736sp'),
	('case2737sop'),
	('case2746wop'),
	('case2746wp'),
	('case2848rte'),
	('case2868rte'),
	('case2869pegase'),
	('case3012wp'),
	('case3120sp'),
	('case3375wp'),
	('case6468rte'),
	('case6470rte'),
	('case6495rte'),
	('case6515rte'),
	('case9241pegase'),
	('case_ACTIVSg2000'),
	('case_ACTIVSg10k');

--CREATE TABLE IF NOT EXISTS linpf.idmap(
--	idint SMALLSERIAL,
--	id VARCHAR(12),
--	PRIMARY KEY(idint)
--);
--
--INSERT INTO linpf.idmap (id) SELECT distinct id from linpf.vt ORDER BY id;
--CREATE UNIQUE INDEX ON linpf.idmap(id);

