.PHONY: docs
docs:
	@julia --project=docs/ --color=yes docs/make.jl
