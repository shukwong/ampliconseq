
version ?= latest
registry = ghcr.io/shukwong

build:
	docker build --tag ${registry}/ampliconseq:${version} .

release: build
	docker push ${registry}/ampliconseq:${version}

singularity: build
	singularity build ampliconseq-${version}.sif docker-daemon://${registry}/ampliconseq:${version}
