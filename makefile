PACKAGE=mcrypt
DOCKER_IMG_NAME=$(PACKAGE)

all: rundocker

builddocker:
	sudo docker build -t $(DOCKER_IMG_NAME) .

rundocker: builddocker
	sudo docker run -it $(DOCKER_IMG_NAME)

