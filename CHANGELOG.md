# Changelog

## [0.11.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.10.2...v0.11.0) (2021-12-21)


### Features

* update to rust-bio-tools 0.33 ([#388](https://www.github.com/IKIM-Essen/uncovar/issues/388)) ([8f66915](https://www.github.com/IKIM-Essen/uncovar/commit/8f669159d7db8dad503cccf505c34a8271795bcd))
* update to varlociraptor 4.9 ([#424](https://www.github.com/IKIM-Essen/uncovar/issues/424)) ([5a7ca5b](https://www.github.com/IKIM-Essen/uncovar/commit/5a7ca5be08afdc6642a0db103c9fec2f6842fc53))

### [0.10.2](https://www.github.com/IKIM-Essen/uncovar/compare/v0.10.1...v0.10.2) (2021-12-19)


### Bug Fixes

* pangolin input for overview table, bamclipper samtools error, include flag for high quality genome report ([#418](https://www.github.com/IKIM-Essen/uncovar/issues/418)) ([9690f12](https://www.github.com/IKIM-Essen/uncovar/commit/9690f12e6725447e035903d3b4eecf31d01b1576))

### [0.10.1](https://www.github.com/IKIM-Essen/uncovar/compare/v0.10.0...v0.10.1) (2021-12-17)


### Bug Fixes

* names in report ([#410](https://www.github.com/IKIM-Essen/uncovar/issues/410)) ([e5dbcbd](https://www.github.com/IKIM-Essen/uncovar/commit/e5dbcbd61fc836e8705a08ac189235c78fd1db65))

## [0.10.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.9.0...v0.10.0) (2021-12-16)


### Features

* complete workflow if no samples pass the quality filter ([#411](https://www.github.com/IKIM-Essen/uncovar/issues/411)) ([55b2dcc](https://www.github.com/IKIM-Essen/uncovar/commit/55b2dcce83c6349b4651c07acfda2cfe885ed0f6))

## [0.9.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.8.0...v0.9.0) (2021-12-15)


### Features

* added variants-over-time-plot to report again ([#405](https://www.github.com/IKIM-Essen/uncovar/issues/405)) ([bd5048a](https://www.github.com/IKIM-Essen/uncovar/commit/bd5048af11d3f1f28a1853a9ec1070fe85ff977b))

## [0.8.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.7.0...v0.8.0) (2021-12-15)


### Features

* add pangolin call table ([#397](https://www.github.com/IKIM-Essen/uncovar/issues/397)) ([953211a](https://www.github.com/IKIM-Essen/uncovar/commit/953211a552f0ce84594e598318cc25a7d406c3cd))
* add rki include/exclude flag ([#319](https://www.github.com/IKIM-Essen/uncovar/issues/319)) ([d2e6d20](https://www.github.com/IKIM-Essen/uncovar/commit/d2e6d20de0429323989d107eac5c4282561b57ce))
* added column in overview table for the clear name of the lineage ([#400](https://www.github.com/IKIM-Essen/uncovar/issues/400)) ([38e34bb](https://www.github.com/IKIM-Essen/uncovar/commit/38e34bba65ab28481466e5897019f060d57596f8))


### Bug Fixes

* change report structure, add masked sequence, pangolin call on polished sequence ([#396](https://www.github.com/IKIM-Essen/uncovar/issues/396)) ([d2ecaaa](https://www.github.com/IKIM-Essen/uncovar/commit/d2ecaaa55efaf916fe380f042758f1560a5b389f))
* rmv pipe from pre-commit ([#401](https://www.github.com/IKIM-Essen/uncovar/issues/401)) ([3222bb2](https://www.github.com/IKIM-Essen/uncovar/commit/3222bb2220ce517ba514f28f416c730ef884d52f))

## [0.7.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.6.0...v0.7.0) (2021-12-14)


### Features

* add pre commit checks ([#398](https://www.github.com/IKIM-Essen/uncovar/issues/398)) ([6982b3d](https://www.github.com/IKIM-Essen/uncovar/commit/6982b3dcaca308840d9a70caec22e79bdb9025fc))


### Bug Fixes

* remove cwd-check in update sample sheet script ([#393](https://www.github.com/IKIM-Essen/uncovar/issues/393)) ([44d49d3](https://www.github.com/IKIM-Essen/uncovar/commit/44d49d34e7ed2ad94002ed9b9e8868ec9665a541))
* trigger docs build on release created ([#391](https://www.github.com/IKIM-Essen/uncovar/issues/391)) ([e240aeb](https://www.github.com/IKIM-Essen/uncovar/commit/e240aeb5f13d38a2ccf14ce1cb2277af197d9117))

## [0.6.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.5.0...v0.6.0) (2021-12-13)


### Features

* update pangolin to 3.1.17 ([#389](https://www.github.com/IKIM-Essen/uncovar/issues/389)) ([23ed39e](https://www.github.com/IKIM-Essen/uncovar/commit/23ed39e6b7471a7f4efc649efb51882cf692e125))
* update varlociraptor to version 4.8 ([#331](https://www.github.com/IKIM-Essen/uncovar/issues/331)) ([0096376](https://www.github.com/IKIM-Essen/uncovar/commit/009637621a3f55b8852c2aa9b428775f1286f40b))


### Bug Fixes

* generation of high quality genomes ([#394](https://www.github.com/IKIM-Essen/uncovar/issues/394)) ([da00975](https://www.github.com/IKIM-Essen/uncovar/commit/da00975893a20144c28cf4185c9263f6bd6c5a7c))
* relative aggregation of samples for variants and lineages over time plots ([#382](https://www.github.com/IKIM-Essen/uncovar/issues/382)) ([ebe556c](https://www.github.com/IKIM-Essen/uncovar/commit/ebe556c3dfa7d1bbf6797d51452d06e53184f679))

## [0.5.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.4.1...v0.5.0) (2021-12-09)


### Features

* **technology:** add ion torrent processing ([#383](https://www.github.com/IKIM-Essen/uncovar/issues/383)) ([288777c](https://www.github.com/IKIM-Essen/uncovar/commit/288777cef21485fb0365c5c552690f2bf94b6e11))

### [0.4.1](https://www.github.com/IKIM-Essen/uncovar/compare/v0.4.0...v0.4.1) (2021-12-07)


### Bug Fixes

* **ci:** adjust commit message of copyright-preamble.yml to fit conventional commits format ([#369](https://www.github.com/IKIM-Essen/uncovar/issues/369)) ([e067895](https://www.github.com/IKIM-Essen/uncovar/commit/e067895adc368414000cba5e4db73edbb3d03b95))
* GISAID lineage extraction [B.1.1.529 (probably).fasta] ([#372](https://www.github.com/IKIM-Essen/uncovar/issues/372)) ([7750862](https://www.github.com/IKIM-Essen/uncovar/commit/7750862a98890937885164b469b5f487ae70488a))

## [0.4.0](https://www.github.com/IKIM-Essen/uncovar/compare/v0.3.1...v0.4.0) (2021-12-02)


### Features

* **technology:** add Oxford Nanopore processing ([#305](https://www.github.com/IKIM-Essen/uncovar/issues/305)) ([55f38f3](https://www.github.com/IKIM-Essen/uncovar/commit/55f38f3f94146355c245f76d6c210dbf351a6eaf))


### Bug Fixes

* add bam index for extract_reads_of_interest ([#366](https://www.github.com/IKIM-Essen/uncovar/issues/366)) ([5f5dd27](https://www.github.com/IKIM-Essen/uncovar/commit/5f5dd27311d1bc27787b2061bba2312a50fbf720))
* changed contigs file for quast-rule to avoid workflow to stop because of low read numbers ([#338](https://www.github.com/IKIM-Essen/uncovar/issues/338)) ([2d0e246](https://www.github.com/IKIM-Essen/uncovar/commit/2d0e2465217470f0e9f089b9c3e240cddfc6336e))
* changed JSON schema validator and updated files [#364](https://www.github.com/IKIM-Essen/uncovar/issues/364) ([75e62f2](https://www.github.com/IKIM-Essen/uncovar/commit/75e62f20b88eaeb19389908e325b264b0f64081d))


### Performance Improvements

* changed chunksize of read json file from gisaid, excluded lineage without name ([#328](https://www.github.com/IKIM-Essen/uncovar/issues/328)) ([b9e5c49](https://www.github.com/IKIM-Essen/uncovar/commit/b9e5c49dbf8910ad0332a73dba0c027bf3545c81))

### [0.3.1](https://www.github.com/IKIM-Essen/uncovar/compare/v0.3.0...v0.3.1) (2021-11-29)


### Bug Fixes

* nested function usage and remove code duplication ([#323](https://www.github.com/IKIM-Essen/uncovar/issues/323)) ([91e3ce9](https://www.github.com/IKIM-Essen/uncovar/commit/91e3ce922132309d3d226a9efceb995ce99b4489))

## [0.3.0](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.5...v0.3.0) (2021-11-05)


### Features

* add aggregation rule for publication plots ([#275](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/275)) ([a93bc18](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/a93bc18bae5976bcdb045d8daab6760e03e3a6ea))
* add more voc mutations to config ([#314](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/314)) ([5098605](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/509860512f819887cfcb15786f127157382d5b29))
* add plot of reads needed for sufficient lineage calling ([#263](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/263)) ([982f997](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/982f997b73fa110babde87cf9bea0406350fb0e7))
* add plot of variants over time ([#279](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/279)) ([66ab121](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/66ab12115814597a569a9662cba51ac679b1a81d))
* add schemes for sample sheet and config file ([#278](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/278)) ([09976b3](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/09976b3b901f224374cc9baa3df1a65eed89822c))
* add standard reference lineages for kallisto based lineage calling ([#280](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/280)) ([6bb32f6](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/6bb32f6cda5cddf6cb598e566bbee9b2af26f190))
* add todo action ([#277](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/277)) ([c1e44c5](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/c1e44c5ee2ad7c9bd2b9cbf5182616a36c50a38d))
* add workflow catalog yml ([#199](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/199)) ([2013533](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/20135333faf2247f2a11e3201d362125de53e3b4))
* comparison of assemblers ([#172](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/172)) ([8233829](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/8233829467b9cfaa45291225d05a8d6e34da044c))


### Bug Fixes

* add "snakemake" in the readme for the Snakemake Workflow Catalog requirements ([#309](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/309)) ([2e5212c](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/2e5212c08604cc4b61704fc07711d00f964078d1))
* changed call config file in update-sample-sheet ([#307](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/307)) ([406aef7](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/406aef7ae411cc13ef9b1df3892922862104d01d))
* shotgun typo in config.yaml  ([#308](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/308)) ([d30f02b](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/d30f02b756699647bfbdcce802a33442bf324a76))
* specify snakemake version used in pangolin env ([#284](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/284)) ([79a12a8](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/79a12a8fbecb0a83f263779880f9384cecf57529))

### [0.2.5](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.4...v0.2.5) (2021-09-23)


### Bug Fixes

* sha tag name ([#273](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/273)) ([b6d1265](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/b6d126582c630d64341ac18f18c982a584022f7e))

### [0.2.4](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.3...v0.2.4) (2021-09-23)


### Bug Fixes

* add log in to quay registry ([#271](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/271)) ([16080de](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/16080de6a68c551698be9be6dbabd441f55befb0))

### [0.2.3](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.2...v0.2.3) (2021-09-23)


### Bug Fixes

* update readme to trigger release ([#269](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/269)) ([19129a9](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/19129a91fae5a068deb55d5243fc1c47953d5d3c))

### [0.2.2](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.1...v0.2.2) (2021-09-22)


### Bug Fixes

* rmv path for buildah ([#265](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/265)) ([2341a61](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/2341a61623d4a905c8ace15363e57db8ada173dd))

### [0.2.1](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.2.0...v0.2.1) (2021-09-22)


### Bug Fixes

* container version tag ([#262](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/262)) ([d4c1ef6](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/d4c1ef6b69439b92932b8398bf1cc8c753f16003))

## [0.2.0](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.8...v0.2.0) (2021-09-22)


### Features

* adding new rule update_sample to workflow ([#254](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/254)) ([3878c1b](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/3878c1b0089396e79dac614c09bf3985d74de96b))


### Bug Fixes

* maximize build space for container image ([#261](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/261)) ([0de2cc9](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/0de2cc92c75bfa4f69ca48b576a57bd688f539b1))

### [0.1.8](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.7...v0.1.8) (2021-09-20)


### Bug Fixes

* add missing path to the Dockerfile job ([#258](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/258)) ([e733d2d](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/e733d2dd1081e632ce39b71cb47d122225011b6d))

### [0.1.7](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.6...v0.1.7) (2021-09-20)


### Bug Fixes

* update snakemake actions for container file ([#256](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/256)) ([cb31e92](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/cb31e92d263dfcdc94675aa3ef8a9274b472654a))

### [0.1.6](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.5...v0.1.6) (2021-09-16)


### Bug Fixes

* change path of docker file ([#251](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/251)) ([1abebcf](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/1abebcff23ef6c4eb82a414a9a9ceadd84cc0c4d))

### [0.1.5](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.4...v0.1.5) (2021-09-15)


### Bug Fixes

* add GISAID var to env ([#248](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/248)) ([e859aed](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/e859aedea76e37594c6a1070827dbc59c3a69b83))

### [0.1.4](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.3...v0.1.4) (2021-09-14)


### Bug Fixes

* change " to ' ([#246](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/246)) ([ffca497](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/ffca4975410e8e1e6de96f058c907e19bc5ed5d5))

### [0.1.3](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.2...v0.1.3) (2021-09-09)


### Bug Fixes

* add extract gisaid env ([#235](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/235)) ([2b94752](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/2b94752704a669d0a29d05eab23bb9e8d7d11fc5))
* improve if statement with release check ([#239](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/239)) ([12cc85f](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/12cc85f672c71361da48b617869582207c5d3f3e))

### [0.1.2](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.1...v0.1.2) (2021-09-06)


### Bug Fixes

* merge push to quay.io into release please ([#231](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/231)) ([e9d782d](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/e9d782dc5db6bd24d097d95f6f945cb02a1ff3aa))
* remove pre von relase action ([#232](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/232)) ([301043f](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/301043f7e1b5a913f3b216eaf66640ae79756f67))

### [0.1.1](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/compare/v0.1.0...v0.1.1) (2021-09-01)


### Bug Fixes

* change trigger of registry release action ([#223](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/223)) ([1464ece](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/1464ece23a689df6f5c9ef0db3a1252b53612cfc))

## 0.1.0 (2021-08-31)


### Bug Fixes

* change release-please-action ([#217](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/217)) ([20af1db](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/20af1dba16d69802a5ed48d453f3a7dc024c26e1))
* change trigger of release please ([#212](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/212)) ([c7a9978](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/c7a9978c8d4439ae77042703b58daa93b416a8a7))
* relase plase indent ([#216](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/216)) ([9140a29](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/9140a2916a34e616c39de4c052567c2dbc51bd3d))
* release command ([#219](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/219)) ([fb84404](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/fb84404bad162ebdc939d8a67ed379ea86d77b75))
* release please command ([#218](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/218)) ([7c64689](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/7c646897d2ff69289b85f24a11bb641628e78c86))
* release please if statement ([#220](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/220)) ([41f41d5](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/41f41d5f21a1f3b1864ee2ad6a1abc8e15a50593))
* release please local PR checkout ([#221](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/221)) ([b348aa7](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/b348aa7ea1eebf00621ee62c52aca655fa8f4c83))
* release please pr tag ([#215](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/215)) ([5d54bb4](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/5d54bb4435d6ac8cdaf696c2b5b991a6ceb49c11))
* release please rmv double dollar ([#222](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/222)) ([0003da9](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/0003da90e8bf65b34ace990576d4e27f5494330c))
* release please step id ([#214](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/214)) ([fa5a9cc](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/fa5a9cccf6930058d169ff2e4a8cbab30bde0ad7))


### Miscellaneous Chores

* release 0.1.0 ([#201](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/issues/201)) ([f211058](https://www.github.com/koesterlab/snakemake-workflow-sars-cov2/commit/f211058ffe0aa8c58b2d35bfb84838ec031b1283))
