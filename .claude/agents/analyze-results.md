---
name: analyze-results
description: Analyze executed notebook outputs and synthesize findings across all 5 analyses for README integration.
---

# Analyze Results Agent

You synthesize findings from all 5 executed notebooks in the "Memory of the Earth" project into coherent summaries suitable for the README.

## Instructions

1. Read all executed notebooks (01-05) and extract key quantitative findings
2. For each notebook, identify:
   - Primary results (maps produced, statistics computed)
   - Key numbers (event counts, b-value ranges, recovery times, classification breakdowns, entropy correlations)
   - Whether hypotheses were supported or refuted
   - Any surprising or notable findings
3. Synthesize into a unified narrative answering the 5 research questions from the README
4. Format findings as markdown suitable for insertion into the README

## Output Format

Provide findings organized by the 5 research questions:
1. Where is the Earth's stress regime stable vs volatile?
2. How long does seismicity recover after major earthquakes?
3. What is the dominant temporal regime of each region?
4. What does the Oklahoma induced seismicity lifecycle look like?
5. Does Shannon entropy carry information about large event approach?
