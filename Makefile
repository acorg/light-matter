.PHONY: check, tcheck, pep8, pyflakes, lint, wc, clean, update-prosite, perf

check: light/_distance.so
	python -m discover -v

tcheck: light/_distance.so
	trial --rterrors test

light/_distance.so: src/distance.py
	python setup.py build_ext -i

pep8:
	find . -path ./light/colors -prune -o \
            -path ./light/performance/data/polymerase/zScores.py -prune -o \
            -path ./light/performance/data/polymerase/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a/zScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_4mtp_a/zScores.py -prune -o \
            -path ./light/performance/data/pdb_4mtp_a/bitScores.py -prune -o \
            -path ./light/performance/data/ha/zScores.py -prune -o \
            -path ./light/performance/data/ha/bitScores.py -prune -o \
            -name '*.py' -print0 | xargs -0 pep8 --ignore=E402
	pep8 --ignore=E201,E241,E501 \
            light/performance/data/ha/zScores.py \
            light/performance/data/polymerase/zScores.py \
            light/performance/data/pdb_2hla_a/zScores.py \
            light/performance/data/pdb_4mtp_a/zScores.py
	pep8 --ignore=E121 \
            light/performance/data/ha/bitScores.py \
            light/performance/data/polymerase/bitScores.py \
            light/performance/data/pdb_2hla_a/bitScores.py \
            light/performance/data/pdb_4mtp_a/bitScores.py

pyflakes:
	find . -path ./light/colors/six.py -prune -o -name '*.py' -print0 | xargs -0 pyflakes 2>&1 | bin/check-pyflakes-output.sh

lint: pep8 pyflakes

wc:
	find . -path ./light/colors -prune -o -path ./test/colors -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name __pycache__ -type d -print0 | xargs -0 rmdir
	find . -name _trial_temp -type d -print0 | xargs -0 rm -r
	rm -f light/*.so
	rm -fr build

perf:
	light/performance/bin/perf.py

performance-data:
	find light/performance/data -name __pycache__ -prune -o -type d -depth 1 -print0 | xargs -0 -n 1 make -C

update-prosite:
	@echo "Downloading Prosite database... this may take a while."
	curl 'ftp://ftp.expasy.org/databases/prosite/prosite.dat' > /tmp/prosite.dat
	bin/prosite-to-json.py /tmp/prosite.dat > data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json
	@echo "New database stored into data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json"
	@echo "You'll need to add the new db to git and update the version number in _DB_FILE in light/landmarks/prosite.py"
	rm /tmp/prosite.dat
