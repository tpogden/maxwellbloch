.DEFAULT_GOAL := all

# Solve -----------------------------------------------------------------------

MBS_JSON = $(wildcard src/mbs*.json)
MBS_QU = $(MBS_JSON:.json=.qu)

src/%.qu: src/%.json
	@echo Solving from $<.
	mbsolve --file $< --recalc

solve: $(MBS_QU)

# Plots -----------------------------------------------------------------------

PLOTS_PY = $(wildcard src/plot_*.py)
PLOTS_PNG = $(PLOTS_PY:.py=.png)

src/plot_%.png: src/plot_%.py src/%.qu
	@echo Plotting $@.
	python $< --png

plots_png: $(PLOTS_PNG)

# Anims -----------------------------------------------------------------------

ANIM_FLAGS = --zoom 1 --fps 30 -c 0.05
ANIMS_PY = $(wildcard src/anim_*.py)
ANIMS_MP4 = $(ANIMS_PY:.py=.mp4)

src/anim_%.mp4: src/anim_%.py src/%.qu
	@echo Animating $@.
	python $< --mp4 $(ANIM_FLAGS)

anims_mp4: $(ANIMS_MP4)

# All -------------------------------------------------------------------------

all: plots_png anims_mp4

# Clean -----------------------------------------------------------------------

clean_qu:
	rm -f src/*.qu

clean_png:
	rm -f src/plot_*.png

clean_mp4:
	rm -f src/anim_*.mp4

clean: clean_qu clean_png clean_mp4
