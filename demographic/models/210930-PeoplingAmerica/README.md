# Infering a four population model that includes the Native American Expansion


We start with the base model [OOA](../210729-ooa/). And add the NAT component.

Copy inferred model and (manually) we added the NAT parameters.

```bash
cp ../210729-ooa/results/ooa_best_fit_model.yml ooa-4pops.yml
```

We add the following parameters to the inferred model:

```
- name: MXB
  description: Native American,  Mexico.
  ancestors: [CHB]
  start_time: 20000
  epochs:
  - {end_time: 10000, start_size: 1000}
  - {end_time: 400, end_size: 50000}
  - {end_time: 0, end_size: 2000}
```
