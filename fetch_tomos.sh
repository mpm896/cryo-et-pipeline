# Matthew Martinez
# Sanofi - US Dept. of Large Molecules Research, Protein Engineering Group, Structural Biology
#
# Fill in values below including your name and the dates of tomograms you want to fetch

DB_DIR=/root/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-data  # Database directory containing datasets. Likely should not change this
NAME="Matthew Martinez"
DATES="2024-07-25"  # Comma separates list with no spaces, in date format YYYY-MM-DD
PIPE_CLI=1	        # 0: Pipe CLI not available
		            # 1. Pipe CLI is available
GRAB_METADATA=1     # 0: Only grab the final tomogram
                    # 1: Grab metadata (i.e. you want to do denoising)


# ============================
# DO NOT MODIFY ANYTHING BELOW 
# ============================

# Get initials
initials=""
for word in $NAME; do 
    initials+=${word:0:1}
done
initials=$(echo $initials | tr '[:upper:]' '[:lower:]')

for date in $(echo $DATES | sed "s/,/ /g"); do 
    # Get DIR prefix 
    DIR=$DB_DIR/$initials$date

    # Find all _rec.mrc files, transfer to dir of same name in the workdir
    if [[ $PIPE_CLI -eq 0 ]]; then
        if [[ $GRAB_METADATA -eq 0 ]]; then
	        find $DIR-* -name "*_rec.mrc" | awk -F "/" '{print "mkdir "$(NF-1)"; rsync --progress -avzhr "$0" "$(NF-1)}'
        elif [[ $GRAB_METADATA -eq 1 ]]; then
            rsync --progress --ignore-existing --exclude={'Frames','Snapshots'} -avzhr $DIR-* .
        fi
    elif [[ $PIPE_CLI -eq 1 ]]; then 
        if [[ $GRAB_METADATA -eq 0 ]]; then
            MAGELLAN_DIR=$(echo $DB_DIR | sed 's/.*\(its\)/\1/g') 
            find $DIR-* -name "*_rec.mrc" \
                | awk -F "/" -v var=$MAGELLAN_DIR '{
                    print "mkdir "$(NF-1)"; pipe storage cp --skip-existing --force cp://"var"/"$(NF-1)"/"$NF" "$(NF-1)"/"
                }' | sh
        elif [[ $GRAB_METADATA -eq 1 ]]; then
            for dataset in $DIR-*; do
                MAGELLAN_DIR=$(echo $dataset | sed 's/.*\(its\)/\1/g')
                pipe storage cp --skip-existing --force -r -e Frames/* -e Snapshots/* cp://$MAGELLAN_DIR $(basename $MAGELLAN_DIR)
            done
        fi
    fi
done