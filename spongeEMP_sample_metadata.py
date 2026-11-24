
# https://github.com/amnona/SpongeEMP/tree/master/sponge_emp/data
file_in = '/Users/songweizhi/Desktop/SpongeEMP/The_sponge_microbiome_project_sample_metadata.txt'

location_set = set()
sponge_associated_sample_num = 0
col_index = dict()
line_num_index = 0
for each_line in open(file_in):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        country = line_split[col_index['country']]
        depth = line_split[col_index['depth']]
        host_common_name = line_split[col_index['host_common_name']]
        host_scientific_name = line_split[col_index['host_scientific_name']]

        # # get the number of sponge_associated_sample
        # if host_common_name not in ['sea water', 'Biofilm', 'marine sediment', 'sea cucumber']:
        #     sponge_associated_sample_num += 1
        #     print('%s\t%s\t%s' % (country, depth, host_scientific_name))

        if country not in ['Palau', 'Netherlands:Netherlands Antilles:Bonaire', 'Klein Curacao', 'The Netherlands', 'Greece', 'Republic of Ireland', 'Norway', 'The Netherlands', 'Curacao', 'United States of America of America', 'Saba Island Region', 'Bonaire Island', 'Curacao Island', 'United Kingdom', 'Gulf of Mexico', 'New Zealand', 'Italy', 'Commonwealth of Puerto Rico', 'Portugal', 'Mexico', 'The Bahamas', 'Canary Islands', 'Israel', 'France', 'Australia', 'Papua New Guinea', 'Tanzania', 'Spain', 'Sweden', 'Bahamas', 'Panama', 'USA', 'Belize', 'Saudi Arabia', 'Croatia']:
            if host_common_name not in ['sea water', 'Biofilm', 'marine sediment', 'sea cucumber']:

                sponge_associated_sample_num += 1
                location_set.add(country)
                print('%s\t%s\t%s' % (country, depth, host_scientific_name))

                # if 'Plakortis simplex' in each_line:
                #     print('%s\t%s\t%s' % (country, depth, host_scientific_name))


print('location_set: %s'                    % location_set)
print('sponge_associated_sample_num: %s'    % sponge_associated_sample_num)


'''
Hyrtios altum					Porifera (Phylum) Demospongiae (Class) Keratosa (Subclass) Dictyoceratida (Order) Thorectidae (Family) Thorectinae (Subfamily) Hyrtios (Genus) Hyrtios altus (Species)
Carteriospongia foliascens		Porifera (Phylum) Demospongiae (Class) Keratosa (Subclass) Dictyoceratida (Order) Thorectidae (Family) Phyllospongiinae (Subfamily) Carteriospongia (Genus) Carteriospongia foliascens (Species)
Pseudoceratina  				Porifera (Phylum) Demospongiae (Class) Verongimorpha (Subclass) Verongiida (Order) Pseudoceratinidae (Family) Pseudoceratina (Genus)
Ianthella basta					Porifera (Phylum) Demospongiae (Class) Verongimorpha (Subclass) Verongiida (Order) Ianthellidae (Family) Ianthella (Genus) Ianthella basta (Species)
Aplysinella    					Porifera (Phylum) Demospongiae (Class) Verongimorpha (Subclass) Verongiida (Order) Aplysinellidae (Family) Aplysinella (Genus)
Agelas  						Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Agelasida (Order) Agelasidae (Family) Agelas (Genus)
Axinella  					    Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Axinellida (Order) Axinellidae (Family) Axinella (Genus)
Rhaphoxya  				        Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Bubarida (Order) Dictyonellidae (Family) Rhaphoxya (Genus)
Acanthella cavernosa			Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Bubarida (Order) Dictyonellidae (Family) Acanthella (Genus) Acanthella cavernosa (Species)
Phakellia fusca					Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Bubarida (Order) Bubaridae (Family) Phakellia (Genus) Phakellia fusca (Species)
Xestospongia     				Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Haplosclerida (Order) Petrosiidae (Family) Xestospongia (Genus)
Callyspongia                    Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Haplosclerida (Order) Callyspongiidae (Family) Callyspongia (Genus)
Stylissa massa					Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Scopalinida (Order) Scopalinidae (Family) Stylissa (Genus) Stylissa massa (Species)
Terpios hoshinota				Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Suberitida (Order) Suberitidae (Family) Terpios (Genus) Terpios hoshinota (Species)
Theonella swinhoei				Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Tetractinellida (Order) Astrophorina (Suborder) Theonellidae (Family) Theonella (Genus) Theonella swinhoei (Species)
Rhabdastrella globostellata		Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Tetractinellida (Order) Astrophorina (Suborder) Ancorinidae (Family) Rhabdastrella (Genus) Rhabdastrella globostellata (Species)
Scleritoderma nodosum			Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Tetractinellida (Order) Spirophorina (Suborder) Scleritodermidae (Family) Scleritoderma (Genus)
Plakortis simplex				Porifera (Phylum) Homoscleromorpha (Class) Homosclerophorida (Order) Plakinidae (Family) Plakortis (Genus) Plakortis simplex (Species)

Aphrocallistes                  Porifera (Phylum) Hexactinellida (Class) Hexasterophora (Subclass) Sceptrulophora (Order) Aphrocallistidae (Family) Aphrocallistes (Genus) Aphrocallistes vastus (Species)
Farrea                          Porifera (Phylum) Hexactinellida (Class) Hexasterophora (Subclass) Sceptrulophora (Order) Farreidae (Family) Farrea (Genus)
Paratimea                       Porifera (Phylum) Demospongiae (Class) Heteroscleromorpha (Subclass) Axinellida (Order) Stelligeridae (Family) Paratimea (Genus) Paratimea alijosensis (Species)           

'''
