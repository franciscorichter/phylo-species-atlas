# site/data

Per-group notes and the few tree files the browser site holds.

```
partitions/<group-slug>/
├── info.yaml                 notes on the group: methods of its reference tree, caveats,
│                             the provenance of its estimated total (text only; the numbers
│                             the page shows come from the paper's files)
├── papers/sources.json       DOI, citation and licence of the papers the notes cite
└── sub-phylos/<slug>/        a tree of a part of the group, not a reference tree
    ├── info.yaml
    └── tree.nwk              identical to standardized/trees/<slug>.nwk where that file
                              exists; otherwise a tree only the site shows (Basidiomycota,
                              Xiphosura) or one derived from the redistributed file
                              (Primates, one tip per species)
```

There is no `partitions/<group-slug>/tree.nwk`: the page reads every reference tree from
`standardized/trees/`. `papers/*.pdf` are not committed (`.gitignore`).
