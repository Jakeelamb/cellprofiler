#!/bin/bash
# Script to mount the easystore removable disk

# Check for USB/external drives
echo "Checking for removable disks..."
lsblk -o NAME,SIZE,TYPE,MOUNTPOINT,FSTYPE,LABEL | grep -E "disk|part"

# Try to find easystore or external USB drives
DEVICES=$(lsblk -o NAME,TYPE,LABEL | grep -E "disk|part" | grep -v "nvme0n1\|zram" | awk '{print $1}')

if [ -z "$DEVICES" ]; then
    echo ""
    echo "No removable disks detected."
    echo "Please:"
    echo "  1. Plug in your USB/external drive"
    echo "  2. Wait a few seconds for it to be detected"
    echo "  3. Run this script again"
    exit 1
fi

echo ""
echo "Found devices: $DEVICES"
echo ""

# Try to auto-mount using udisks2 (if available)
if command -v udisksctl &> /dev/null; then
    for dev in $DEVICES; do
        if [[ $dev =~ ^sd[a-z]$ ]] || [[ $dev =~ ^sd[a-z][0-9]+$ ]]; then
            echo "Attempting to mount /dev/$dev..."
            udisksctl mount -b "/dev/$dev" 2>&1
            MOUNT_POINT=$(lsblk -o MOUNTPOINT -n "/dev/$dev" | head -1)
            if [ -n "$MOUNT_POINT" ] && [ "$MOUNT_POINT" != "" ]; then
                echo "Successfully mounted at: $MOUNT_POINT"
                echo ""
                echo "You can now access your files at:"
                echo "  $MOUNT_POINT"
                exit 0
            fi
        fi
    done
fi

# Fallback: manual mount instructions
echo "To mount manually, find your device and run:"
echo "  sudo mkdir -p /mnt/easystore"
echo "  sudo mount /dev/sdX1 /mnt/easystore"
echo ""
echo "Replace sdX1 with your actual device (check with 'lsblk')"

