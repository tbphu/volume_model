VENV=virtualenv
VENV_DIR=${PWD}/env
PYTHON=${VENV_DIR}/bin/python
PIP=${VENV_DIR}/bin/pip

default: virtualenv install

virtualenv:
	${VENV} --no-site-packages -p python2 ${VENV_DIR};

install: virtualenv 	
	${PIP} install -r requirements.txt
