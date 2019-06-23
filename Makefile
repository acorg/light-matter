.PHONY: check, tcheck, pyflakes, lint, wc, clean, update-prosite, perf

# If you are on OS X, you'll need the GNU version of find. If you're using
# brew, run brew install findutils (which installs gfind in /usr/local/bin).
FIND := $(shell which gfind || which find)

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)


check: light/_distance.so
	python -m discover -v

tcheck: light/_distance.so
	trial --rterrors test

light/_distance.so: src/distance.py
	python setup.py build_ext -i

pycodestyle:
	$(FIND) . -path ./light/colors -prune -o \
            -path ./light/performance/data/polymerase/zScores.py -prune -o \
            -path ./light/performance/data/polymerase/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a/zScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_4mtp_a/zScores.py -prune -o \
            -path ./light/performance/data/pdb_4mtp_a/bitScores.py -prune -o \
            -path ./light/performance/data/ha/zScores.py -prune -o \
            -path ./light/performance/data/ha/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_4ph0_a/zScores.py -prune -o \
            -path ./light/performance/data/pdb_4ph0_a/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a_against_polymerase/bitScores.py -prune -o \
            -path ./light/performance/data/pdb_2hla_a_against_ha/bitScores.py -prune -o \
            -path ./test/test_bin_score.py -prune -o \
            -path ./test/test_binScoreTemplate.py -prune -o \
            -name '*.py' -print0 | $(XARGS) -0 pycodestyle --ignore=E402
	pycodestyle --ignore=E201,E241,E501 \
            light/performance/data/ha/zScores.py \
            light/performance/data/polymerase/zScores.py \
            light/performance/data/pdb_2hla_a/zScores.py \
            light/performance/data/pdb_4mtp_a/zScores.py \
            light/performance/data/pdb_4ph0_a/zScores.py
	pycodestyle --ignore=E121 \
            light/performance/data/ha/bitScores.py \
            light/performance/data/polymerase/bitScores.py \
            light/performance/data/pdb_2hla_a/bitScores.py \
            light/performance/data/pdb_4mtp_a/bitScores.py \
            light/performance/data/pdb_4ph0_a/bitScores.py \
            light/performance/data/pdb_2hla_a_against_polymerase/bitScores.py \
            light/performance/data/pdb_2hla_a_against_ha/bitScores.py
	pycodestyle --ignore=E501 \
            test/test_bin_score.py \
            test/test_binScoreTemplate.py

pyflakes:
	$(FIND) . -path ./light/colors/six.py -prune -o -name '*.py' -print0 | $(XARGS) -0 pyflakes 2>&1 | bin/check-pyflakes-output.sh

lint: pycodestyle

wc:
	$(FIND) . -path ./light/colors -prune -o -path ./test/colors -prune -o -name '*.py' -print0 | $(XARGS) -0 wc -l

clean:
	$(FIND) . \( -name '*.pyc' -o -name '*~' \) -print0 | $(XARGS) -0 rm
	$(FIND) . -name __pycache__ -type d -print0 | $(XARGS) -0 rmdir
	$(FIND) . -name _trial_temp -type d -print0 | $(XARGS) -0 rm -r
	rm -f light/*.so
	rm -fr build

perf:
	light/performance/bin/perf.py

performance-data:
	$(FIND) light/performance/data -maxdepth 1 -mindepth 1 \( \! -name __pycache__ \) -type d -print0 | $(XARGS) -0 -n 1 make -C

update-prosite:
	@echo "Downloading Prosite database... this may take a while."
	curl 'ftp://ftp.expasy.org/databases/prosite/prosite.dat' > /tmp/prosite.dat
	bin/prosite-to-json.py /tmp/prosite.dat > data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json
	@echo "New database stored into data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json"
	@echo "You'll need to add the new db to git and update the version number in _DB_FILE in light/landmarks/prosite.py"
	rm /tmp/prosite.dat
