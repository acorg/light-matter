.PHONY: check, tcheck, pep8, pyflakes, lint, wc, clean, update-prosite

check:
	python -m discover -v

tcheck:
	trial --rterrors test

pep8:
	find . -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find . -path './light/colors/six.py' -prune -o -name '*.py' -print0 | xargs -0 pyflakes 2>&1 | bin/check-pyflakes-output.sh

lint: pep8 pyflakes

wc:
	find . -path './light/colors' -prune -o -path './test/colors' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
	find . -name '_trial_temp' -type d -print0 | xargs -0 rm -r

update-prosite:
	@echo "Downloading Prosite database... this may take a while."
	curl 'ftp://ftp.expasy.org/databases/prosite/prosite.dat' > /tmp/prosite.dat
	bin/prosite-to-json.py /tmp/prosite.dat > data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json
	@echo "New database stored into data/prosite-`egrep -m 1 '^CC   Release ' /tmp/prosite.dat | awk '{print $$3}'`.json"
	@echo "You'll need to add the new db to git and update the version number in _DB_FILE in light/landmarks/prosite.py"
	rm /tmp/prosite.dat
