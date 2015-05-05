.PHONY: check, tcheck, pep8, pyflakes, lint, wc, clean

check:
	python -m discover -v

pep8:
	find . -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find . -path './light/colors/six.py' -prune -o -name '*.py' -print0 | xargs -0 pyflakes

lint: pep8 pyflakes

wc:
	find . -path './light/colors' -prune -o -path './test/colors' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
